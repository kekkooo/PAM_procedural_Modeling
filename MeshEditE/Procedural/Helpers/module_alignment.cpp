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
void build_manifold_kdtree( Manifold& m, vector< VertexID > &selected, kd_tree &tree)
{
    for( VertexID v : selected )
    {
        tree.insert( m.pos(v), v);
    }
    tree.build();    
}

/// builds a pseudo-random transformation matrix
void generate_random_transform( HMesh::Manifold &host, HMesh::Manifold &module, CGLA::Mat4x4d &t )
{
    // calculate bounding sphere of the two manifolds
    Vec3d  host_centroid, module_centroid;
    double host_radius, module_radius;
    Procedural::Geometry::bsphere( module, module_centroid,   module_radius );
    Procedural::Geometry::bsphere( host,   host_centroid,     host_radius   );
    
    // choose a random edge as translation direction and another as rotation axis
    size_t  tr_edge_number  = CGLA::gel_rand() % host.no_halfedges();
    size_t  rot_edge_number = CGLA::gel_rand() % module.no_halfedges();
    int     count           = 0;
    HalfEdgeID tr_candidate;
    for( auto it = host.halfedges_begin(); it != host.halfedges_end() && count < tr_edge_number; ++it, ++count )
    {
        tr_candidate = *it;
    }
    assert( tr_candidate != InvalidHalfEdgeID );
    count = 0;
    HalfEdgeID rot_candidate;
    for( auto it = module.halfedges_begin(); it != module.halfedges_end() && count < rot_edge_number; ++it, ++count )
    {
        rot_candidate = *it;
    }
    assert( rot_candidate != InvalidHalfEdgeID );
    
    cout << "translataion -> " << tr_candidate << " # rotation " << rot_candidate << endl;
    
    // get the actual vectors
    Vec3d tr_dir    = host.pos( host.walker( tr_candidate ).vertex( )) -
    host.pos( host.walker( tr_candidate ).opp().vertex( ));
    Vec3d rot_axis  = module.pos( module.walker( rot_candidate ).vertex( )) -
    module.pos( module.walker( rot_candidate ).opp().vertex( ));
    tr_dir += rot_axis;
    tr_dir.normalize();
    tr_dir *= ( host_radius );
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

/// applu pseudo-transformation to the poles of the module
void transform_module_poles( HMesh::Manifold &host, HMesh::Manifold &module, std::map< HMesh::VertexID, CGLA::Vec3d > &new_pos)
{
    CGLA::Mat4x4d t;
    generate_random_transform( host, module, t );
    
    for( auto pole : module.vertices( ))
    {
        if( is_pole( module, pole )){
            new_pos[pole] = t.mul_3D_point( module.pos( pole ));
        }
    }
}

/// find correspondances between module and
void match_module_to_host( Manifold &host, Manifold &module, kd_tree &tree,
                           vector< VertexID > &host_p, vector< VertexID > &module_p )
{
    for( auto vid : module.vertices( ))
    {
        if( is_pole( module, vid )){
            double      distance = numeric_limits<double>::max();
            Vec3d       foundVec;
            VertexID    foundID;
            // here there could be problems because I should be sure that there will not be
            // two poles assigned to the same candidate vertex on the host
            assert( tree.closest_point( module.pos( vid ), distance, foundVec, foundID ));
            
            module_p.push_back( vid );
            host_p.push_back( foundID );
        }
    }
    assert(module_p.size() == host_p.size());
}

EdgeCost get_glueings( Manifold &host, Manifold &module, vector<VertexID> &host_vertices, vector<VertexID> &module_poles,
                       size_t no_glueings, vector<Match> &matches )
{
    EdgeCost cost = get_best_poles_subset(host, host_vertices, module, module_poles, matches, no_glueings);
    return cost;
}
            

void apply_optimal_alignment( Manifold &host, Manifold &module, vector<Match> &matches )
{
    Mat4x4d R, T;
    vector< VertexID > host_v, module_p;
    for( auto m : matches )
    {
        module_p.push_back(m.first);
        host_v.push_back(m.second);
    }
    svd_rigid_motion( module, module_p, host, host_v, R, T );
    
    Mat4x4d t = T * R;
//    cout << "rotation "     << R << endl;
//    cout << "translation "  << T << endl;
//    cout << t << endl;
    for( auto mv : module.vertices() ) { module.pos( mv ) = t.mul_3D_point( module.pos( mv )); }
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

void add_necessary_poles( Manifold &host, vector<VertexID> &selected,
                          Manifold &module, vector<VertexID> &poles, vector<VertexID> &new_ids )
{
    HMesh::VertexAttributeVector<int> vertex_selection;
    // for each selected vertex
    for( int i = 0; i < selected.size(); ++i )
    {
        if( !is_pole( host, selected[i] ))
        {
            // clean the selection
            for( auto v : host.vertices( )) { vertex_selection[v] = 0; }
            // calculate the right size accordingly to the valence of the corresponding pole
            /// TODOOOOOOO
                                        size_t size = -1;
            assert(size != -1);
            // call add branch
            VertexID poleID = add_branch( host, selected[i], size, vertex_selection );
            assert( poleID != InvalidVertexID );
            // save the returned value into new_ids
            new_ids.push_back( poleID );
        }
        else
        {
            new_ids.push_back(selected[i]);
        }
    }
}


            
void AddModule( Manifold &host, Manifold &module, size_t no_glueings )
{    
    vector<matches_and_cost> matches_vector;

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
}}}