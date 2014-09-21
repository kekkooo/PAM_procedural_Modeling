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
#include <algorithm>
#include <queue>
#include <sstream>
#include <GEL/CGLA/Mat4x4d.h>
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Operations/geometric_operations.h>
#include <MeshEditE/test.h>
//#include <polarize.h>

using namespace HMesh;
using namespace CGLA;
using namespace std;
using namespace Geometry;


namespace Procedural{
    namespace Matching{
        

        void console_call( Manifold& me_active_mesh )
        {
            string path = "/Users/francescousai/Documents/Dottorato/Visiting/Results/tests/";
            
            for( int i = 0; i < 100; ++i)
            {
                Manifold module;
                build( me_active_mesh );
                build( module );
                double angle = match( me_active_mesh, module );
                
                stringstream oss;
                oss << path;

                if ( angle < M_PI_2 ) oss << "discard/";

                oss << "tri" << i <<  "_" << angle << ".obj";
                obj_save( oss.str(), me_active_mesh );
            }
        }
        
        // loads the choosen step, loads the module and glues it to the active mesh
        void console_call2( Manifold& me_active_mesh, int baseID )
        {
            // initialize the meshes
            Manifold module;
            build( me_active_mesh, baseID );
//            build( me_active_mesh, baseID );
            build( module );
//            polar_subdivide( me_active_mesh, 1 );
            
            // scale active
            double scaling_factor = 0.75 + ((double)( gel_rand() % 100 )/200.0);
            Mat4x4d scale = scaling_Mat4x4d(Vec3d( scaling_factor, scaling_factor, scaling_factor ));
            for( auto v : me_active_mesh.vertices())
            {
                me_active_mesh.pos(v) = scale.mul_3D_point(me_active_mesh.pos(v));
            }
//            return;
            // calculate distances
            DistanceVector combined_dist;
            HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( me_active_mesh );
            Procedural::Structure::LabelJunctions( me_active_mesh, edge_info );
            Procedural::Geometry::distance_from_poles_and_junctions ( me_active_mesh, edge_info, combined_dist, false   );
            
            // select from the me_active_mesh points that have the right combined distance
            vector<VertexID> selected;
            select_candidate_vs( me_active_mesh, selected, combined_dist );

            // move the module in a random position in space
            transform( module, me_active_mesh );
            set< VertexID > module_poles, fresh_module_ids;
            // add the module manifold to me_active_mesh ( the manifold contains 2 separate components
            // save the vertexID of the module's poles and all of its vertex ids in the new manifold
            copy( module, me_active_mesh, module_poles, fresh_module_ids );
            
            // choose the best matching pole
            VertexID    candidate_v = InvalidVertexID, candidate_pole = InvalidVertexID;
            double      min_dist    = find_best_matching( me_active_mesh, selected, module_poles,
                                                          candidate_v, candidate_pole );
            
            assert( candidate_pole != InvalidVertexID );
            assert( candidate_v    != InvalidVertexID );
            
            // remove candidates from usable poles and usable vertices
            module_poles.erase( candidate_pole );
            selected.erase( std::find(selected.begin(), selected.end(), candidate_v));

            // take the orignial vertex normal
            Vec3d   vn  = Geometry::vertex_normal( me_active_mesh, candidate_v );
            candidate_v = add_pole_if_necessary( me_active_mesh, candidate_v, 3 );
            
            // save intermediate result
            string path = "/Users/francescousai/Documents/Dottorato/Visiting/Results/tests/";
            stringstream oss;
            oss << path << "before.obj";
            obj_save(oss.str(), me_active_mesh );
            
            // ALIGNMENT
            // take candidate_pole and candidate_v normals
            align_module( me_active_mesh, module, candidate_pole, vn, fresh_module_ids );
            
            // save intermediate result
            oss = stringstream();
            oss << path << "intermediate.obj";
            obj_save(oss.str(), me_active_mesh );

            Vec3d translation_dir = me_active_mesh.pos( candidate_v ) - me_active_mesh.pos( candidate_pole );
            Mat4x4d tr_pole_to_v = translation_Mat4x4d( translation_dir );
            Mat4x4d push         = translation_Mat4x4d( vn );
            Mat4x4d t2           = push * tr_pole_to_v;
            for( auto mv : fresh_module_ids )
            {
                me_active_mesh.pos(mv) = t2.mul_3D_point(me_active_mesh.pos(mv));
            }
            
            // now that the first pole is fixed, we need to find another glueing point
            // find the rotation that brings the other pole to that vertex
            // track the movement of the vertex to that pole
            Vec3d old_pole_pos = me_active_mesh.pos(candidate_pole);
            Procedural::Operations::Structural::glue_poles( me_active_mesh, candidate_v, candidate_pole );
            
            // save intermediate result
            oss = stringstream();
            oss << path << "intermediate1.obj";
            obj_save(oss.str(), me_active_mesh );


            VertexID    other_pole      = InvalidVertexID,
                        other_candidate = InvalidVertexID;
            
            find_best_matching( me_active_mesh, selected, module_poles, other_candidate, other_pole );
            Vec3d other_candidate_pos = me_active_mesh.pos(other_candidate);
            Vec3d other_vn  = Geometry::vertex_normal( me_active_mesh, other_candidate );
            other_candidate = add_pole_if_necessary( me_active_mesh, other_candidate, 3 );
//            Procedural::Operations::Geometric::extrude_pole( me_active_mesh, other_candidate, other_candidate_pos - me_active_mesh.pos(other_candidate), false, 1.0 );
            
//            align_module( me_active_mesh, module, other_pole, other_vn, fresh_module_ids );
            
            Vec3d   pn              = Geometry::vertex_normal( me_active_mesh, other_pole );
            Vec3d   rotation_axis   = CGLA::cross( other_vn, pn );
            double  rotation_angle  = std::acos( CGLA::dot( other_vn, pn ));
            double  pace            = rotation_angle / 2.0;
            double  r;
            Vec3d   pivot;
            bsphere_of_selected( me_active_mesh, fresh_module_ids, pivot, r );
            Mat4x4d tr_origin       = translation_Mat4x4d( -pivot );
            Mat4x4d rot             = get_rotation_mat4d( rotation_axis, pace );
            Mat4x4d tr_back         = translation_Mat4x4d( pivot );
            Mat4x4d t               = tr_back * rot * tr_origin;
            
            cout << "angle : " << rotation_angle << " # pace: " << pace << endl;
            
            Vec3d control_point = t.mul_3D_point( me_active_mesh.pos(other_pole ));
            vector< Vec3d > points;
            bezier( me_active_mesh.pos(other_pole), control_point, other_candidate_pos, 7, points );
            
            
//            Procedural::Operations::Geometric::add_ring_around_pole( me_active_mesh, other_pole, 1.0 );


            for( Vec3d cbp : points )
            {
                Vec3d current = me_active_mesh.pos( other_pole );
                double radius = Procedural::Geometry::ring_mean_radius(me_active_mesh, me_active_mesh.walker(other_pole).next().halfedge());
                Procedural::Operations::Geometric::extrude_pole( me_active_mesh, other_pole, cbp - current, true, 1.0 );
                double new_radius = Procedural::Geometry::ring_mean_radius(me_active_mesh, me_active_mesh.walker(other_pole).next().halfedge());
                Procedural::Operations::Geometric::scale_ring_radius(me_active_mesh, me_active_mesh.walker(other_pole).next().halfedge(), radius/ new_radius, true );
            }
            
            assert( is_pole( me_active_mesh, other_pole ));

            
            
            
            // save intermediate result
            oss = stringstream();
            oss << path << "intermediate2.obj";
            obj_save(oss.str(), me_active_mesh );

            assert( is_pole( me_active_mesh, other_candidate ));
            Procedural::Operations::Structural::glue_poles( me_active_mesh, other_candidate, other_pole );

            
            
            
            // find a rotation around the normal of  candidate_v  that alings in a good way the vertices
            // choose a vertex on the ring connected to candidate_v
//            VertexID cvvr       = me_active_mesh.walker(candidate_v).next().vertex();
//            Vec3d    ccvr_pos   = me_active_mesh.pos(cvvr);
//            Vec3d    cvn        = Procedural::Geometry::vertex_normal(me_active_mesh, candidate_v);
//            Walker   w          = me_active_mesh.walker( candidate_pole );
//            
//            double min_angle = numeric_limits<double>::max();
//            while( !w.full_circle())
//            {
//                double angle =
//                    std::acos( CGLA::dot(cvn, me_active_mesh.pos(w.vertex()) - ccvr_pos));
//                if( angle < min_angle ) min_angle = angle;
//                w = w.circulate_vertex_ccw();
//                         
//            }
            
            
            // build and apply the rotation matrix with axis = vcn and angle = min_angle/2
//            Vec3d p_normal = Procedural::Geometry::vertex_normal(me_active_mesh, candidate_pole);
//            rot = get_rotation_mat4d( p_normal, -min_angle/2.0);
//            bsphere_of_selected( me_active_mesh, fresh_ids, centroid, radius );
//            tr_origin = translation_Mat4x4d( centroid * ( -1.0 ));
//            tr_back   = translation_Mat4x4d( centroid );
//            t =  tr_back * rot * tr_origin;
//            for( auto mv : fresh_ids )
//            {
//                me_active_mesh.pos(mv) = t.mul_3D_point(me_active_mesh.pos(mv));
//            }
            
//            Procedural::Operations::Structural::glue_poles( me_active_mesh, candidate_v, candidate_pole );
//            alt_glue_poles( me_active_mesh, candidate_v, candidate_pole );
            
            oss = stringstream();
            oss.clear();
            oss << path << "after.obj";
            obj_save(oss.str(), me_active_mesh );
        }
        
        double match( Manifold& destination, Manifold& module )
        {
            double angle = 0.0;
            
            transform( module, destination );
            kd_tree tree;
            build_mesh_kdtree( destination, tree );
            
            set< VertexID > module_poles;
            set< VertexID > _;
            copy( module, destination, module_poles, _ );
            cout << " there are " << module_poles.size() << " poles " << endl;
            
            destination.cleanup();
            
            VertexID pole_to_glue;
            VertexID candidate = find_nearest_match( destination, module_poles, tree, pole_to_glue );
            // calculate angle between the two vertex normals
            Vec3d pole_to_glue_normal   = Procedural::Geometry::vertex_normal( destination, pole_to_glue );
            Vec3d candidate_normal      = Procedural::Geometry::vertex_normal( destination, candidate );
            // calculate the angle
            angle = std::acos( CGLA::dot(pole_to_glue_normal, candidate_normal ));
            cout << " angle : " << angle << endl;
            if ( angle < M_PI_2 ) return angle;
            
            
            
            assert( candidate != InvalidVertexID );
            if( !is_pole( destination, candidate ))
            {
                // find all the poles and after check which is the new one.
                HMesh::VertexAttributeVector<int> ring;
                for( auto v : destination.vertices()) { ring[v] = 0; }
                ring[candidate] = 1;
                set<VertexID> old_poles, new_poles;
                pole_set( destination, old_poles );
                Procedural::Operations::Structural::add_branch( destination, candidate, 3, ring );
                pole_set( destination, new_poles );
                for( auto pole : new_poles )
                {
                    if( old_poles.count( pole ) == 0 )
                        candidate = pole;
                }
            }
            Procedural::Operations::Structural::glue_poles( destination, candidate, pole_to_glue );

            return angle;
        }
        
        
        void copy( Manifold& source, Manifold& destination,
                   set< VertexID > &source_poles_in_dest, set< VertexID > &fresh_vertex_IDs )
        {
            // save old pole VertexID
            set< VertexID > old_poles;
            pole_set( destination, old_poles );
            // save old vertex IDs
            set< VertexID > old_dest_IDS;
            for( auto v : destination.vertices() ) old_dest_IDS.insert( v );
            
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
                if( is_pole( destination, vid ))
                {
                    if( old_poles.count(vid) == 0 )
                        source_poles_in_dest.insert( vid );
                }
            }
                        
            for( auto vid : destination.vertices( ))
            {
                if( old_dest_IDS.count( vid ) == 0 )
                    fresh_vertex_IDs.insert( vid );
            }

            
//            set_difference( current_IDs.begin(), current_IDs.end(), old_dest_IDS.begin(), old_dest_IDS.end(), fresh_vertex_IDs.begin() );
        }
        
        void pole_set ( Manifold& m, set< VertexID > &poles )
        {
            for( auto vid : m.vertices() )
            {
                if( is_pole( m, vid ))
                {
                    poles.insert( vid );
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
        
        
        void build( Manifold& m, int iter )
        {
            m.clear();
            stringstream oss;
            oss << "/Users/francescousai/Documents/Dottorato/Visiting/Results/sep09/ex7/test_match";
            if ( iter > 0 ) oss << iter;
            oss << ".obj";
            obj_load(oss.str(), m );
        }

        
        void bsphere( Manifold& m, Vec3d& centroid, double& radius )
        {
            Vec3d pmin, pmax;
            if(m.no_vertices()==0)
                return;
            VertexIDIterator v = m.vertices_begin();
            pmin = pmax = m.pos(*v);
            ++v;
            for(; v != m.vertices_end(); ++v){
                pmin = v_min(m.pos(*v), pmin);
                pmax = v_max(m.pos(*v), pmax);
            }
            
            Manifold::Vec rad = ( pmax - pmin ) * 0.5f;
            centroid = pmin + rad;
            radius = rad.length();
        }
        
        
        void bsphere_of_selected( Manifold& m, set<VertexID> selected, Vec3d& centroid, double& radius )
        {
            Vec3d pmin, pmax;
            if(m.no_vertices()==0)
                return;
            VertexIDIterator v = m.vertices_begin();
            pmin = pmax = m.pos(*v);
            ++v;
//            for(; v != m.vertices_end(); ++v){
            for( auto v : selected )
            {
                pmin = v_min(m.pos(v), pmin);
                pmax = v_max(m.pos(v), pmax);
            }
            
            Manifold::Vec rad = ( pmax - pmin ) * 0.5f;
            centroid = pmin + rad;
            radius = rad.length();
        }
        
        
        void transform ( Manifold& source, Manifold &destination )
        {
            // calculate bounding sphere of the two manifolds
            Vec3d  source_centroid, dest_centroid;
            double source_radius, dest_radius;
            bsphere( source,        source_centroid,    source_radius );
            bsphere( destination,   dest_centroid,      dest_radius   );
            
            // choose a random edge as translation direction and another as rotation axis
            size_t tr_edge_number = CGLA::gel_rand() % destination.no_halfedges();
            size_t rot_edge_number = CGLA::gel_rand() % source.no_halfedges();
            int count = 0;
            HalfEdgeID tr_candidate;
            for( auto it = destination.halfedges_begin(); it != destination.halfedges_end() && count < tr_edge_number; ++it, ++count )
            {
                tr_candidate = *it;
            }
            assert( tr_candidate != InvalidHalfEdgeID );
            count = 0;
            HalfEdgeID rot_candidate;
            for( auto it = source.halfedges_begin(); it != source.halfedges_end() && count < rot_edge_number; ++it, ++count )
            {
                rot_candidate = *it;
            }
            assert( rot_candidate != InvalidHalfEdgeID );
            
            cout << "translataion -> " << tr_candidate << " # rotation " << rot_candidate << endl;
            
            // get the actual vectors
            Vec3d tr_dir    = destination.pos( destination.walker( tr_candidate ).vertex( ))
                            - destination.pos( destination.walker( tr_candidate ).opp().vertex( ));
            Vec3d rot_axis  = source.pos( source.walker( rot_candidate ).vertex( ))
                            - source.pos( source.walker( rot_candidate ).opp().vertex( ));
            tr_dir += rot_axis;
            tr_dir.normalize();
            tr_dir *= ( dest_radius );
            double  rot_x       = ( gel_rand() % 70 ) / 100.0,
                    rot_y       = ( gel_rand() % 70 ) / 100.0,
                    rot_z       = ( gel_rand() % 70 ) / 100.0,
                    tr_scale    = (( gel_rand() % 25 ) / 100 ) + 1.0;

            auto rot =
                CGLA::rotation_Mat4x4d( XAXIS, rot_x  ) *
                CGLA::rotation_Mat4x4d( YAXIS, rot_y  ) *
                CGLA::rotation_Mat4x4d( ZAXIS, rot_z  );
            cout << rot << endl;
            
            Vec3d translation = tr_dir * tr_scale;
            CGLA::Mat4x4d t = rot * CGLA::translation_Mat4x4d( translation );
        
            for( auto vid : source.vertices( ))
            {
                source.pos(vid) =   t.mul_3D_point( source.pos( vid ));
            }
        }
        
        // this should also discard points that are too distant from the module
        void select_candidate_vs ( Manifold& m, vector< VertexID >& vs, DistanceVector& dist )
        {
            for( VertexID v : m.vertices())
            {
                // magic number!
                if( dist[v] > 3 || is_pole(m, v) ) vs.push_back( v );
            }

        }
        
        double find_best_matching  ( Manifold& m, vector< VertexID > valid_vs,
                                     set< VertexID > module_poles,
                                     VertexID &candidate_v, VertexID &candidate_pole )
        {
            double      min_dist        = numeric_limits<double>::max();
            candidate_v     = InvalidVertexID;
            candidate_pole  = InvalidVertexID;
            for( VertexID v : valid_vs )
            {
                for( VertexID pole : module_poles )
                {
                    double dist = ( m.pos( v ) - m.pos( pole ) ).length();
                    if( dist < min_dist )
                    {
                        min_dist        = dist;
                        candidate_v     = v;
                        candidate_pole  = pole;
                    }
                }
            }
            return min_dist;
        }
        
        
        void align_module ( Manifold& m, Manifold& module, VertexID& pole, Vec3d axis_to_align_with,
                            set< VertexID > points_to_move)
        {
            Vec3d   pn              = Geometry::vertex_normal( m, pole );
            Vec3d   rotation_axis   = CGLA::cross( axis_to_align_with, pn );
            double  rotation_angle  = std::acos( CGLA::dot( axis_to_align_with, pn ));
            Mat4x4d rot             = get_rotation_mat4d( rotation_axis, rotation_angle);
            // find center of the module mesh
            Vec3d  centroid;
            double radius;
            bsphere( module, centroid, radius );
            Mat4x4d tr_origin = translation_Mat4x4d( -centroid );
            Mat4x4d tr_back   = translation_Mat4x4d( centroid );
            Mat4x4d t =  tr_back * rot * tr_origin;
            
            cout << rot         << endl << endl;
            cout << tr_origin   << endl << endl;
            cout << tr_back     << endl << endl;
            cout << t           << endl << endl;
            
            assert( points_to_move.size() > 0 );
            for( auto mv : points_to_move )
            {
                m.pos( mv ) = t.mul_3D_point( m.pos( mv ));
            }
        }
        
        VertexID add_pole_if_necessary( Manifold& m, VertexID candidate_v, int size )
        {
            // if the selected vertex is not a pole
            if( !is_pole( m, candidate_v ))
            {
                // add the new pole in the selected vertex location
                HMesh::VertexAttributeVector<int> ring;
                for( auto v : m.vertices()) { ring[v] = 0; }
                ring[candidate_v] = 1;
                
                set<VertexID> old_poles, new_poles;
                pole_set( m, old_poles );
                Procedural::Operations::Structural::add_branch( m, candidate_v, 3, ring );
                pole_set( m, new_poles );
                
                // find all the poles and after check which is the new one.
                for( auto pole : new_poles )
                {
                    if( old_poles.count( pole ) == 0 ) { candidate_v = pole; }
                }
                assert( is_pole( m, candidate_v ));
            }
            Procedural::Operations::Geometric::flatten_pole( m, candidate_v );
            return  candidate_v;

        }
    }
}
