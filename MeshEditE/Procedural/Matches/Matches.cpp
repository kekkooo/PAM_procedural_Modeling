//
//  Matches.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 04/09/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "Matches.h"
#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <polarize.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/obj_load.h>
#include <set>
#include <queue>
#include <sstream>
#include <GEL/CGLA/Mat4x4d.h>
#include <MeshEditE/Procedural/Operations/structural_operations.h>

using namespace HMesh;
using namespace CGLA;
using namespace std;
using namespace Geometry;

namespace Procedural{
    namespace Matching{
        

        
        void console_call( Manifold& me_active_mesh )
        {
            Manifold module;
            build( me_active_mesh );
            build( module );
            transform(module);
            kd_tree tree;
            build_mesh_kdtree( me_active_mesh, tree );
            
            set< VertexID > module_poles;
            copy( module, me_active_mesh, module_poles );
            cout << " there are " << module_poles.size() << " poles " << endl;

            me_active_mesh.cleanup();
            
            VertexID pole_to_glue;
            VertexID candidate = find_nearest_match( me_active_mesh, module_poles, tree, pole_to_glue );
            assert( candidate != InvalidVertexID );
            if( !is_pole( me_active_mesh, candidate ))
            {
                HMesh::VertexAttributeVector<int> ring;
                for( auto v : me_active_mesh.vertices()) { ring[v] = 0; }
                ring[candidate] = 1;
                Procedural::Operations::Structural::add_branch( me_active_mesh, candidate, 3, ring );
            }
            Procedural::Operations::Structural::glue_poles(me_active_mesh, candidate, pole_to_glue );
            
        }
        
        double match( Manifold& destination, Manifold& module )
        {
            double angle = 0.0;
            
            
            
            
            
            return angle;
        }
        
        
        void copy( Manifold& source, Manifold& destination, set< VertexID > &source_poles_in_dest )
        {
            // save old pole VertexID
            set< VertexID > old_poles;
            for( auto vid : destination.vertices() )
            {
                if( is_pole(destination, vid )) old_poles.insert( vid );
            }
            
            map< FaceID, vector< FaceID > > face_steps_map;
            map< FaceID, FaceID >    old_to_new_faces;

            for( auto f : source.faces( ))
            {
                vector< Vec3d > vs;
                Walker w = source.walker(f);
                while(!w.full_circle())
                {
                    if( face_steps_map.count(f) == 0 )
                        face_steps_map[f] = vector< FaceID >();
                    face_steps_map[f].push_back( w.opp().face() );
                    // prendo il prev, perché altrimenti parto una rotazione avanti
                    vs.push_back( source.pos( w.prev().vertex( )));
                    w = w.next();
                }
                FaceID new_f = destination.add_face( vs );
                old_to_new_faces[f] = new_f;
            }
            
            for( auto f : source.faces( ))
            {
                FaceID new_f   = old_to_new_faces[f];
                auto neighbors = face_steps_map[f];
                Walker fw = destination.walker( new_f );
                for( int steps = 0; steps < neighbors.size(); ++steps )
                {
//                    if( fw.opp().face() != InvalidFaceID ) continue;
                    // find the neighor face
                    FaceID old_f_neighbor = neighbors[steps];
                    FaceID new_neighbor   = old_to_new_faces[old_f_neighbor];
                    // find N that is the index at face_steps_map[old_f_neighbor that gives f
                    int i = 0;
                    for(  ; i < neighbors.size() && face_steps_map[old_f_neighbor][i] != f; ++i );
                    assert( face_steps_map[old_f_neighbor][i] == f );
                    Walker neighbor_w = destination.walker(new_neighbor);
                    for( int j = 0; j < i; ++j ) { neighbor_w = neighbor_w.next(); }
                    destination.stitch_boundary_edges(neighbor_w.opp().halfedge(), fw.opp().halfedge());

                    fw = fw.next();
                }
            }
            
            
            // add mesh to the other
            
            destination.cleanup();
            // find the difference in poles
            cout << " there are " << destination.no_vertices() << " verts " << endl;
            for( auto vid : destination.vertices() )
            {
                if( is_pole(destination, vid ))
                {
                    if( old_poles.count(vid) == 0 )
                        source_poles_in_dest.insert( vid );
                }
            }
        }
        
        
        VertexID find_nearest_match( Manifold& source, set< VertexID > &poles, kd_tree &tree, VertexID &choosen_pole )
        {
            VertexID candidate      = InvalidVertexID;
            double   min_distance   = numeric_limits<double>::max();
            for( VertexID pole : poles )
            {
                double      distance = numeric_limits<double>::max();
                Vec3d       foundVec;
                VertexID    foundID;
                
                bool have_found = tree.closest_point( source.pos(pole), distance, foundVec, foundID );
                if( have_found )
                {
                    if( distance < min_distance )
                    {
                        min_distance = distance;
                        candidate    = foundID;
                        choosen_pole = pole;
                    }
                }
            }
            return candidate;
        }
        
        
        void build_mesh_kdtree( Manifold& m, kd_tree &tree )
        {
            for( VertexID v : m.vertices() )
            {
                tree.insert(m.pos(v), v);
            }
            tree.build();
        }
        
        
        void build( Manifold& m )
        {
            stringstream oss;
            oss << "/Users/francescousai/Documents/Dottorato/Visiting/Results/test_match.obj";
            obj_load(oss.str(), m );
        }
        
        void transform ( Manifold& m )
        {
            Vec3d translation( 10, 10, 10);
            CGLA::Mat4x4d t = CGLA::translation_Mat4x4d( translation );
        
            for( auto vid : m.vertices( ))
            {
                m.pos(vid) = t.mul_3D_point( m.pos( vid ));
            }
            
        }
    }
}
