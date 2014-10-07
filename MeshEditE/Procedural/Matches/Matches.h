//
//  Matches.h
//  MeshEditE
//
//  Created by Francesco Usai on 04/09/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Matches__
#define __MeshEditE__Matches__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <polarize.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/obj_load.h>
#include <set>
#include <sstream>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
//#include <MeshEditE/Procedural/Helpers/svd_rigid_motion.h>
#include <MeshEditE/Procedural/Helpers/svd_alignment.h>

using namespace HMesh;
using namespace CGLA;
using namespace std;
using namespace Geometry;

namespace Procedural{
    namespace Matching{
        
        typedef VertexAttributeVector< Procedural::Geometry::DistanceMetrics >  DistanceVector;
        typedef KDTree< Vec3d, VertexID > kd_tree;
        
        void        console_call        ( Manifold& me_active_mesh );
        void        console_call2       ( Manifold& me_active_mesh, int baseID = 0 );
        void        align_to_selected   ( Manifold& m, vector< VertexID > selected );
        
        void        copy                ( Manifold& source, Manifold& destination,
                                          set< VertexID > &source_poles_in_dest,
                                          set< VertexID > &fresh_vertex_IDs );
        void        copy                ( Manifold& source, Manifold& destination,
                                          set< VertexID > &source_poles_in_dest,
                                          set< VertexID > &fresh_vertex_IDs,
                                          map< VertexID, VertexID > source_to_dest_poles_ID);

        double      match               ( Manifold& destination, Manifold& module );
        VertexID    find_nearest_match  ( Manifold& source, set< VertexID > &poles, kd_tree &tree, VertexID &choosen_pole );
        void        build_mesh_kdtree   ( Manifold& m, kd_tree &tree );
        void        build_mesh_kdtree   ( Manifold& m, set< VertexID > &selected, kd_tree &tree);
        void        build               ( Manifold& m, int iter = 0 );
        void        transform           ( Manifold& source, Manifold& destination );
        void        bsphere             ( Manifold& m, Vec3d& centroid, double& radius );
        void        bsphere_of_selected ( Manifold& m, set<VertexID> selected, Vec3d& centroid, double& radius );
        void        pole_set            ( Manifold& m, set< VertexID > &poles );
        void        select_candidate_vs ( Manifold& m, vector< VertexID > &vs, DistanceVector& dist );
        double      find_best_matching  ( Manifold& m, vector< VertexID > valid_vs,
                                          set< VertexID > module_poles,
                                          VertexID &candidate_v, VertexID &candidate_pole );
        void        align_module        ( Manifold& m, Manifold& module, VertexID& pole,
                                          Vec3d axis_to_align_with, set< VertexID > points_to_move );
        VertexID    add_pole_if_necessary( Manifold& m, VertexID candidate_v, int size );

        
        
            
    }}

    #endif /* defined(__MeshEditE__Matches__) */
