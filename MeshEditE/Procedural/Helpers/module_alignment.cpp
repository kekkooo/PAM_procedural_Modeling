//
//  module_alignment.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 02/11/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "module_alignment.h"
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Matches/graph_match.h>
#include <polarize.h>

using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Geometry;
using namespace Procedural::Geometry;
using namespace Procedural::GraphMatch;

namespace Procedural {
    namespace Helpers{
        namespace ModuleAlignment{
            
void build_manifold_kdtree( Manifold& m, vector< VertexID > &selected, kd_tree &tree)
{
    for( VertexID v : selected )
    {
        tree.insert( m.pos(v), v);
    }
    tree.build();    
}
            
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

void transform_module_poles( HMesh::Manifold &host, HMesh::Manifold &module, std::map< HMesh::VertexID, CGLA::Vec3d > &new_pos)
{
    CGLA::Mat4x4d t;
    generate_random_transform( host, module, t );
    
    for( auto vid : module.vertices( ))
    {
        if( is_pole( module, vid )){
            new_pos[vid] = t.mul_3D_point( module.pos( vid ));
        }
    }
}
            
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
            host_p.push_back(   foundID );
        }
    }
}

void get_glueings( Manifold &host, Manifold &module, vertex_match &pole_to_host_vertex,
                   int no_glueings, vector< VertexID > selected_poles )
{
    
    
}



            
void AddModule( Manifold &host, Manifold &module, size_t no_glueings )
{

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