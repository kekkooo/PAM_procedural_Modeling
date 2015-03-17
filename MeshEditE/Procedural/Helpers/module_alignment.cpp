//
//  module_alignment.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 02/11/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "module_alignment.h"

#include <unordered_set>

#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <MeshEditE/Procedural/Helpers/manifold_copy.h>
#include <polarize.h>
#include "Test.h"
#include <random>


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

/// builds a pseudo-random transformation Tr that brings the module
/// in a possibly collision-free position outside the host's bounding sphere
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
    HalfEdgeID rot_candidate = *m.halfedges_begin();
    for( auto it = m.halfedges_begin(); it != m.halfedges_end() && count < rot_edge_number; ++it, ++count )
    {
        rot_candidate = *it;
    }
    assert( rot_candidate != InvalidHalfEdgeID );
    
    cout << "translataion -> " << tr_candidate << " # rotation " << rot_candidate << endl;
    
    // get the actual vectors
    Vec3d tr_dir    = vec_from_edge( m, tr_candidate );
    Vec3d rot_axis  = vec_from_edge( m, rot_candidate );
    
    tr_dir    += rot_axis;
//    tr_dir.normalize();
//          tr_dir    *= ( host_radius );

    // What the hell I'm doing here???
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
    // try to avoid intersections
    Vec3d transformed_centroid  = t.mul_3D_point( module_centroid );
    // this does not work if the centroid is the same e.g. both on the origin
    Vec3d centroid_dir          = transformed_centroid - host_centroid;
//    put the module outside the hosts bsphere
    if( centroid_dir.length() < host_radius + module_radius )
    {
        double length = host_radius + module_radius - centroid_dir.length();
        centroid_dir.normalize();
        centroid_dir *= length;
        t = CGLA::translation_Mat4x4d( centroid_dir ) * t;
    }

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
    set<VertexID> assigned_vertices;
    for( auto id_and_pos : module_poles_positions)
    {
        assert( is_pole( m, id_and_pos.first ));
        double      distance    = numeric_limits<double>::max();
        Vec3d       foundVec;
        VertexID    foundID = InvalidVertexID, second_choiceID = InvalidVertexID;
        // here there could be problems because I should be sure that there will not be
        // two poles assigned to the same candidate vertex on the host
        // closest point should always return true
        bool        have_found  = tree.closest_point( id_and_pos.second , distance, foundVec, foundID );
        cout << id_and_pos.first << " # " << id_and_pos.second << "# matched : " << foundID << " # dist : " << distance << endl;
        assert( have_found );
        assert( foundID != InvalidVertexID );
        if( assigned_vertices.count( foundID ) == 0 )
        {
            assigned_vertices.insert( foundID );
            pole_to_host_vertex[id_and_pos.first] = foundID;
        }
        else
        {
            find_second_closest( m, tree, foundID, second_choiceID, assigned_vertices );
            assert( second_choiceID != InvalidVertexID );
            assigned_vertices.insert( second_choiceID );
            pole_to_host_vertex[id_and_pos.first] = second_choiceID;
        }
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
//    return;
    save_intermediate_result(m, TEST_PATH , 1);
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
void add_necessary_poles( Manifold &m, vector<Match> &matches, set<VertexID> &host_vs, set<VertexID> &module_vs )
{
    vector<VertexID>                    host_poles_to_remap, host_poles_remapped;
    set<VertexID>                       module_vs_remapped,  host_vs_remapped;
    HMesh::VertexAttributeVector<int>   vertex_selection;
    IDRemap                             remap;
    // aggiungo i poli e salvo il loro id
    for( int i = 0; i < matches.size(); ++i )
    {
        if( !is_pole( m, matches[i].second ))
        {
            // clean the selection
            for( auto v : m.vertices( )) { vertex_selection[v] = 0; }
            // calculate the right size accordingly to the valence of the corresponding pole
            /// TODOOOOOOO
#pragma message("size here is arbitrary")
            int size = 2;
//            int size = -1;
//            assert( size != -1 );
            VertexID poleID = add_branch( m, matches[i].second, size, vertex_selection );
            assert( poleID != InvalidVertexID );
            matches[i].second = poleID;
//            host_poles_to_remap.push_back( poleID );
        }
//        else
//        {
//            host_poles_to_remap.push_back( matches[i].second );
//        }
    }
    m.cleanup( remap );
    // re-align the saved IDs
//    for( VertexID pole : host_poles_to_remap )
//    {
//        assert( remap.vmap.count(pole) > 0 );
//        host_poles_remapped.push_back( remap.vmap[pole] );
//    }
    for( auto p : remap.vmap )
    {
        if( p.second == InvalidVertexID ) { continue; }
        if( module_vs.count( p.first ) > 0 )    { module_vs_remapped.insert( p.second ); }
        else                                    { host_vs_remapped.insert( p.second );   }
    }
    
    // change
    for( int i = 0; i < matches.size(); ++i )
    {
        matches[i].first  = remap.vmap[matches[i].first];
        matches[i].second = remap.vmap[matches[i].second];
    }
//    selected    = std::move( host_poles_remapped );
    module_vs   = std::move( module_vs_remapped );
    host_vs     = std::move( host_vs_remapped );
    
}


            
void AddModule( Manifold &host, Manifold &module, size_t no_glueings, vector<Match> &result )
{
    // initialize random generator
    static std::mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31,
    0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253> rrrr;
    rrrr.seed( time(0) ); rrrr(); rrrr(); rrrr(); rrrr(); rrrr(); rrrr(); rrrr();
    CGLA::gel_srand( rrrr() );

    vector<matches_and_cost>    matches_vector;
    kd_tree                     tree;
    set<VertexID>               host_IDs, module_IDs;
    vector<match_info>          proposed_matches;
    
    add_manifold( host, module, host_IDs, module_IDs );
    build_manifold_kdtree( host, host_IDs, tree );
    // for now it will be just for one time
    std::cout << "while filling " << endl;
    for( int i = 0; i < 10; ++i )
    {
        map<VertexID, Vec3d>    transformed_module_poles;
        vertex_match            module_to_host;
        vector<Match>           current_matches, best_matches;
        match_info              mi;
        
        mi.random_transform =
            transform_module_poles( host, host_IDs, module_IDs, transformed_module_poles );
        
//        for( auto p : transformed_module_poles )
//        {
//            cout << "########################" << endl;
//            cout << "id  : " << p.first << endl;
//            cout << "bef : " << host.pos(p.first) << endl;
//            cout << "aft : " << p.second<< endl;
//            cout << "########################" << endl;
//        }
//#define save_and_continue
#ifdef save_and_continue
        assert(false);
        for( VertexID vid : module_IDs )
        {
            host.pos( vid ) = mi.random_transform.mul_3D_point( host.pos( vid ));
        }
        save_intermediate_result( host, TEST_PATH, 10 + i );
//        continue;
#endif
        
        
        match_module_to_host( host, tree, transformed_module_poles, module_to_host );
//        continue;
        for( auto pole_and_vertex : module_to_host )
        {
//            cout << pole_and_vertex.first << ", " << pole_and_vertex.second << endl;
            current_matches.push_back( make_pair( pole_and_vertex.first, pole_and_vertex.second ));
        }
        EdgeCost c = get_best_subset( host, current_matches, best_matches, no_glueings );

        graph_print( cout, c ) << endl;;
        mi.cost     = c;
        mi.matches  = best_matches;
        proposed_matches.push_back( mi );
        graph_print( cout, proposed_matches[i].cost ) << endl;;
    }
    std::cout << "after filling " << endl;
    assert( proposed_matches.size() > 0 );
    // find the best solution between the proposed ones
#pragma message( "all the matches have the same cost - I guess it is related to the RANDOM transform")
    size_t      selected = 0;
    EdgeCost    max_cost = proposed_matches[0].cost;
    graph_print( cout, max_cost ) << endl;
    for( int i = 1; i < proposed_matches.size(); ++i )
    {
        graph_print( cout, proposed_matches[i].cost ) << endl;;
        if( proposed_matches[i].cost < max_cost )
        {
            selected = i;
            max_cost = proposed_matches[i].cost;
        }
    }
    // FROM NOW ON, THE OPERATIONS WILL USE THE SELECTED MATCHES
    // remember that there is a return in the middle of code
    apply_optimal_alignment( host, module_IDs, proposed_matches[selected] );
    save_intermediate_result(host, TEST_PATH, 2);
    align_module_normals_to_host( host, module_IDs, proposed_matches[selected].matches );
    save_intermediate_result(host, TEST_PATH, 3);
    add_necessary_poles( host, proposed_matches[selected].matches, host_IDs, module_IDs );
    glue_matches( host, proposed_matches[selected].matches );
    

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
    result = std::move( proposed_matches[selected].matches );
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
        if( assigned.count(w.vertex())) { continue; }   // skip all assigned vertices
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
        assert( it != ids_and_dists.end( ));
        done    = ( assigned.count( it->first ) == 0 );
    }
    second_closest = it->first;
    assigned.insert( it->first );    
}
            
            
void align_module_normals_to_host( Manifold &m, set<VertexID> &module_IDs, vector<Match> &matches )
{
    CGLA::Vec3d host_vec( 0 ), module_vec( 0 ), centroid( 0 );
    double      _;
    for( Match match : matches )
    {
        centroid += m.pos( match.first );
        Vec3d mn = vertex_normal( m, match.first );
        Vec3d hn = vertex_normal( m, match.second );
        mn.normalize();
        hn.normalize();
        module_vec  += mn;
        host_vec    -= hn;
    }
    module_vec.normalize();
    host_vec.normalize();
    centroid /= matches.size();
    // need to carefully choose which centroid I should use.
//    bsphere( m, module_IDs, centroid, _ );
    Mat4x4d t = get_alignment_for_2_vectors( module_vec, host_vec, centroid );
    for( VertexID v : module_IDs )
    {
        m.pos(v) = t.mul_3D_point( m.pos( v ));
    }
}
            
void glue_matches( HMesh::Manifold &m, std::vector<Procedural::GraphMatch::Match> &matches )
{
    for( Match match : matches )
    {
        Procedural::Operations::Structural::glue_poles( m, match.first, match.second );
    }
    
}

/*=========================================================================*
 *  R A N D _ R O T A T I O N      Author: Jim Arvo, 1991                  *
 *  Modified by                            Francesco Usai, 2015            *
 *                                                                         *
 *  This routine maps three values (x[0], x[1], x[2]) in the range [0,1]   *
 *  into a 3x3 rotation matrix, M.  Uniformly distributed random variables *
 *  x0, x1, and x2 create uniformly distributed random rotation matrices.  *
 *  To create small uniformly distributed "perturbations", supply          *
 *  samples in the following ranges                                        *
 *                                                                         *
 *      x[0] in [ 0, d ]                                                   *
 *      x[1] in [ 0, 1 ]                                                   *
 *      x[2] in [ 0, d ]                                                   *
 *                                                                         *
 * where 0 < d < 1 controls the size of the perturbation.  Any of the      *
 * random variables may be stratified (or "jittered") for a slightly more  *
 * even distribution.                                                      *
 *                                                                         *
 *=========================================================================*/
Mat4x4d random_rotation_matrix_arvo( float x1, float x2, float x3 )
{
    Mat4x4d M; // = rotation_Mat4x4d(XAXIS, M_PI * 2.0 * x_1 );
    
    float theta = x1 * M_PI * 2.0; /* Rotation about the pole (Z).      */
    float phi   = x2 * M_PI * 2.0; /* For direction of pole deflection. */
    float z     = x3 * 2.0;        /* For magnitude of pole deflection. */

    /* Compute a vector V used for distributing points over the sphere  */
    /* via the reflection I - V Transpose(V).  This formulation of V    */
    /* will guarantee that if x[1] and x[2] are uniformly distributed,  */
    /* the reflected points will be uniform on the sphere.  Note that V */
    /* has length sqrt(2) to eliminate the 2 in the Householder matrix. */
    
    float r  = sqrt( z );
    float Vx = sin( phi ) * r;
    float Vy = cos( phi ) * r;
    float Vz = sqrt( 2.0 - z );

    /* Compute the row vector S = Transpose(V) * R, where R is a simple */
    /* rotation by theta about the z-axis.  No need to compute Sz since */
    /* it's just Vz.                                                    */
    
    float st = sin( theta );
    float ct = cos( theta );
    float Sx = Vx * ct - Vy * st;
    float Sy = Vx * st + Vy * ct;
    
    /* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
    /* is equivalent to V S - R.                                        */
    
    float m_00 = Vx * Sx - ct,
          m_01 = Vx * Sy - st,
          m_02 = Vx * Vz,
          m_03 = 0;
    Vec4d r0 = Vec4d( m_00, m_01, m_02, m_03 );
    
    float m_10 = Vy * Sx + st,
          m_11 = Vy * Sy - ct,
          m_12 = Vy * Vz,
          m_13 = 0;
    Vec4d r1 = Vec4d( m_10, m_11, m_12, m_13 );
    
    float m_20 = Vz * Sx,
          m_21 = Vz * Sy,
          m_22 = 1.0 - z,   /* This equals Vz * Vz - 1.0 */
          m_23 = 0;
    Vec4d r2 = Vec4d( m_20, m_21, m_22, m_23 );
    Vec4d r3 = Vec4d( 0.0, 0.0, 0.0, 1.0 );
    
    M = Mat4x4d( r0, r1, r2, r3 );
    return M;
}
  
}}}