//
//  module_alignment.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 02/11/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "module_alignment.h"
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <polarize.h>
#include <unordered_set>
#include <MeshEditE/Procedural/Helpers/manifold_copy.h>

using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Geometry;
using namespace Procedural::Geometry;
using namespace Procedural::GraphMatch;
using namespace Procedural::Operations::Structural;

namespace Procedural {
    namespace Helpers{
        namespace ModuleAlignment{

/// prepare the spatial index
void build_manifold_kdtree( Manifold& m, set< VertexID > &selected, kd_tree &tree)
{
    for( VertexID v : selected )
    {
        tree.insert( m.pos(v), v);
    }
    tree.build();    
}

/// builds a pseudo-random transformation matrix
void generate_random_transform( HMesh::Manifold &m, const set<VertexID> &host_vs,
                                const set<VertexID> &module_vs, CGLA::Mat4x4d &t )
{
    // calculate bounding sphere of the two manifolds
    Vec3d  host_centroid, module_centroid;
    double host_radius, module_radius;
    Procedural::Geometry::bsphere( m, module_vs, module_centroid,   module_radius );
    Procedural::Geometry::bsphere( m, host_vs,   host_centroid,     host_radius   );
    
    // choose a random edge as translation direction and another as rotation axis
    size_t  tr_edge_number  = CGLA::gel_rand() % m.no_halfedges();
    size_t  rot_edge_number = CGLA::gel_rand() % m.no_vertices();
    int     count           = 0;
    HalfEdgeID tr_candidate;
    for( auto it = m.halfedges_begin(); it != m.halfedges_end() && count < tr_edge_number; ++it, ++count )
    {
        tr_candidate = *it;
    }
    assert( tr_candidate != InvalidHalfEdgeID );
    count = 0;
    HalfEdgeID rot_candidate;
    for( auto it = m.halfedges_begin(); it != m.halfedges_end() && count < rot_edge_number; ++it, ++count )
    {
        rot_candidate = *it;
    }
    assert( rot_candidate != InvalidHalfEdgeID );
    
    cout << "translataion -> " << tr_candidate << " # rotation " << rot_candidate << endl;
    
    // get the actual vectors
    Vec3d tr_dir    = m.pos( m.walker( tr_candidate ).vertex( )) -
                      m.pos( m.walker( tr_candidate ).opp().vertex( ));
    Vec3d rot_axis  = m.pos( m.walker( rot_candidate ).vertex( )) -
                      m.pos( m.walker( rot_candidate ).opp().vertex( ));
          tr_dir    += rot_axis;
    tr_dir.normalize();
          tr_dir    *= ( host_radius );

    double  rot_x       = ( gel_rand() % 70 ) / 100.0,
            rot_y       = ( gel_rand() % 70 ) / 100.0,
            rot_z       = ( gel_rand() % 70 ) / 100.0,
            tr_scale    = (( gel_rand() % 25 ) / 100 ) + 1.0;
    
    auto rot = CGLA::rotation_Mat4x4d( XAXIS, rot_x  ) *
               CGLA::rotation_Mat4x4d( YAXIS, rot_y  ) *
               CGLA::rotation_Mat4x4d( ZAXIS, rot_z  );
    cout << rot << endl;
    
    Vec3d translation = tr_dir * tr_scale;
    t = rot * CGLA::translation_Mat4x4d( translation );
}

/// m contains both the host and the module
Mat4x4d transform_module_poles( HMesh::Manifold &m, const set<VertexID> &host_vs, const set<VertexID> &module_vs,
                                std::map< HMesh::VertexID, CGLA::Vec3d > &new_pos)
{
    CGLA::Mat4x4d t;
    generate_random_transform( m, host_vs, module_vs, t );
    for( auto pole : module_vs )
    {
        if( is_pole( m, pole )){
            new_pos[pole] = t.mul_3D_point( m.pos( pole ));
        }
    }
    return t;
}

/// find correspondances between module and host
void match_module_to_host( Manifold &m, kd_tree &tree, map<VertexID, Vec3d> &module_poles_positions,
                           vertex_match &pole_to_host_vertex )
{
    for( auto id_and_pos : module_poles_positions)
    {
        assert( is_pole( m, id_and_pos.first ));
        double      distance = numeric_limits<double>::max();
        Vec3d       foundVec;
        VertexID    foundID;
        // here there could be problems because I should be sure that there will not be
        // two poles assigned to the same candidate vertex on the host
        // closest point should always return true
        bool        have_found = tree.closest_point( id_and_pos.second , distance, foundVec, foundID );
        assert( have_found );
        
        pole_to_host_vertex[id_and_pos.first] = foundID;
    }
}
            

void apply_optimal_alignment( Manifold &m, const set<VertexID> &module_vs, match_info &choosen_match_info )
{
    Mat4x4d R, T;
    vector< VertexID > host_v, module_p;
    for( auto m : choosen_match_info.matches )
    {
        module_p.push_back(m.first);
        host_v.push_back(m.second);
    }
    
    // before calculating the optimal rigid motion, transform the module using the given initial transformation
    // that is the random transformation used.
    for( auto mv : module_vs ) { m.pos( mv ) = choosen_match_info.random_transform.mul_3D_point(m.pos( mv )); }
    svd_rigid_motion( m, module_p, m, host_v, R, T );
    Mat4x4d t = T * R;
//    cout << "rotation "     << R << endl;
//    cout << "translation "  << T << endl;
//    cout << t << endl;
    for( auto mv : module_vs ) { m.pos( mv ) = t.mul_3D_point( m.pos( mv )); }
}
            
void best_configuration( vector<matches_and_cost> matches_and_costs, matches_and_cost &choosen_cost )
{
    size_t      index       = 0;
    EdgeCost    best_cost   = matches_and_costs[0].second;
    for( int i = 1; i < matches_and_costs.size(); ++i)
    {
        if( matches_and_costs[i].second < best_cost )
        {
            best_cost   = matches_and_costs[i].second;
            index       = i;
        }
    }
    choosen_cost = matches_and_costs[index];
}

/// If some of the selected vertices from the host are not poles, then they must be converted into poles.
/// Since doing this changes the structure of the mesh, we need to clean up the Manifold
/// and update the sets of vertices' IDs of the host and the module
void add_necessary_poles( Manifold &m, vector<VertexID> &selected, set<VertexID> &host_vs, set<VertexID> &module_vs )
{
    vector<VertexID>                    host_poles_to_remap, host_poles_remapped;
    set<VertexID>                       module_vs_remapped,  host_vs_remapped;
    HMesh::VertexAttributeVector<int>   vertex_selection;
    IDRemap                             remap;
    // aggiungo i poli e salvo il loro id
    for( int i = 0; i < selected.size(); ++i )
    {
        if( is_pole( m, selected[i] ))
        {
            // clean the selection
            for( auto v : m.vertices( )) { vertex_selection[v] = 0; }
            // calculate the right size accordingly to the valence of the corresponding pole
            /// TODOOOOOOO
            int size = -1;
            assert( size != -1 );
            VertexID poleID = add_branch( m, selected[i], size, vertex_selection );
            assert( poleID != InvalidVertexID );
            host_poles_to_remap.push_back( poleID );
        }
        else
        {
            host_poles_to_remap.push_back( selected[i] );
        }
    }
    m.cleanup( remap );
    // re-align the saved IDs
    for( VertexID pole : host_poles_to_remap )
    {
        assert( remap.vmap.count(pole) > 0 );
        host_poles_remapped.push_back( remap.vmap[pole] );
    }
    for( auto p : remap.vmap )
    {
        if( p.second == InvalidVertexID ) { continue; }
        if( module_vs.count( p.first ) > 0 )    { module_vs_remapped.insert( p.second ); }
        else                                    { host_vs_remapped.insert( p.second );   }
    }
    
    selected    = std::move( host_poles_remapped );
    module_vs   = std::move( module_vs_remapped );
    host_vs     = std::move( host_vs_remapped );
    
}


            
void AddModule( Manifold &host, Manifold &module, size_t no_glueings )
{    
    vector<matches_and_cost>    matches_vector;
    kd_tree                     tree;
    set<VertexID>               host_IDs, module_IDs;
    vector<match_info>          proposed_matches;
    
    add_manifold( host, module, host_IDs, module_IDs );
    build_manifold_kdtree( host, host_IDs, tree );
    // for now it will be just for one time
    {
        map<VertexID, Vec3d>    transformed_module_poles;
        vertex_match            module_to_host;
        vector<Match>           current_matches, best_matches;
        match_info              mi;
        
        mi.random_transform =
            transform_module_poles( host, host_IDs, module_IDs, transformed_module_poles );
        match_module_to_host( host, tree, transformed_module_poles, module_to_host );
        for( auto pole_and_vertex : module_to_host )
        {
            current_matches.push_back( make_pair( pole_and_vertex.first, pole_and_vertex.second ));
        }
        EdgeCost c = get_best_subset( host, current_matches, best_matches, no_glueings );

        mi.cost     = c;
        mi.matches  = best_matches;
        proposed_matches.push_back( mi );
    }
    assert(proposed_matches.size() > 0);
    // find the best solution between the proposed ones
    size_t      selected = 0;
    EdgeCost    max_cost = proposed_matches[0].cost;
    for( int i = 1; i < proposed_matches.size(); ++i )
    {
        if( proposed_matches[i].cost < max_cost )
        {
            selected = i;
            max_cost = proposed_matches[i].cost;
        }
    }
    // FROM NOW ON, THE OPERATIONS WILL USE THE SELECTED MATCHES
    apply_optimal_alignment( host, module_IDs, proposed_matches[selected] );
    

    /* How this should work
     -) build a kdtree of the host
     -) for N times
        -)  compute a random transformation of the module (that brings it outside the host)
            it only apply the transformation to its poles
        -)  match each pole of the module to its closest point on the host
        -)  use the complete graph utilities in order to get the subset ( with given cardinality )
            that is the best fit between module and host
        -)  collect the configuration and its COST
     -) choose the configuration with smallest cost
     -) use svd rigid motion utility to align the module to the host given the matching between
        the choosen poles and vertices
     -) add a pole to the host in those choosen vertices that are not poles
     -) glue each choosen module's pole to its correspondent host's pole.
     -) smooth the added skeleton ( this is something that needs a little more work )
     */
}
            
bool ID_and_dist_comparer( ID_and_dist &l, ID_and_dist &r ) { return ( l.second < r.second ); }
            
void find_second_closest( const Manifold &m, const kd_tree &tree, const VertexID &closest,
                          VertexID &second_closest, set<VertexID> &assigned)
{
    assert( assigned.count(closest) > 0 );
    // consider closest
    Vec3d       cp          = m.pos( closest );
    // find the nearest of the vertices in its 1-ring
    double      min_dist    = numeric_limits<double>::max();
    for( Walker w = m.walker( closest ); !w.full_circle(); w = w.circulate_vertex_ccw())
    {
        double dist = ( cp - m.pos( w.vertex( ))).length( );
        if( dist < min_dist ){ min_dist = dist; }
    }
    // take that distance as radius and use it with tree.in_sphere with center in closest
    double              radius          = min_dist * 1.1;
    vector<VertexID>    in_sphere_ID;
    vector<Vec3d>       in_sphere_points;
    IDs_and_dists       ids_and_dists;

    tree.in_sphere( cp, radius, in_sphere_points, in_sphere_ID );
    // take the nearest point
    assert( in_sphere_ID.size() == in_sphere_points.size());
    for( int i = 0; i < in_sphere_points.size(); ++i)
    {
        ids_and_dists.push_back( make_pair(in_sphere_ID[i], (cp - in_sphere_points[i]).length()));
    }
    // sort in ascending order accordingly to the distance
    std::sort( ids_and_dists.begin(), ids_and_dists.end(), ID_and_dist_comparer );
    IDs_and_dists::iterator it      = ids_and_dists.begin();
    bool                    done    = ( assigned.count( it->first ) == 0 );
    while ( !done )
    {
        ++it;
        done    = ( assigned.count( it->first ) == 0 );
    }
    second_closest = it->first;
    assigned.insert( it->first );    
}

}}}