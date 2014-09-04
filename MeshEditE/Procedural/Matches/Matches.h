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

using namespace HMesh;
using namespace CGLA;
using namespace std;
using namespace Geometry;

namespace Procedural{
    namespace Matching{
        
        typedef KDTree< Vec3d, VertexID > kd_tree;
        
        void        copy                ( Manifold& source, Manifold& destination, set< VertexID > &source_poles_in_dest );
        double      match               ( Manifold& destination, Manifold& module );
        VertexID    find_nearest_match  ( Manifold& source, set< VertexID > &poles, kd_tree &tree, VertexID &choosen_pole );
        void        build_mesh_kdtree   ( Manifold& m, kd_tree &tree );
        void        build               ( Manifold& m );
        void        transform           ( Manifold& m );
        void        console_call        ( Manifold& me_active_mesh );
        
}}

#endif /* defined(__MeshEditE__Matches__) */
