//
//  test_random_rotation.h
//  MeshEditE
//
//  Created by Francesco Usai on 11/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__test_random_rotation__
#define __MeshEditE__test_random_rotation__

#include <stdio.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Mat4x4d.h>
#include <vector>
#include <map>
#include <random>
#include <set>

namespace Tests{
    typedef std::vector<HMesh::VertexID>                    VertexList;
    typedef std::map<HMesh::VertexID, CGLA::Vec3d>          VertexPosMap;
    typedef std::set<HMesh::VertexID>                       VertexSet;

    class RandomRotationTest{
        
        public :
        
        static  RandomRotationTest&  getInstance( );
                void                 setHost( HMesh::Manifold &manifold );
                void                 setModule( HMesh::Manifold &manifold );
                void                 releaseMesh( );
                void                 rotate( );
                void                 translateOutsideBSphere( );
                void                 resetOriginalModulePosition();
                
        private :
        RandomRotationTest();
        // C++ 11
        // =======
        // We can use the better technique of deleting the methods
        // we don't want.
                RandomRotationTest(RandomRotationTest const&)   = delete;
        void    operator = (RandomRotationTest const&)          = delete;
        
        
        // reference to mesh
        HMesh::Manifold *m;
        
        VertexSet           H_vertices,
                            M_vertices;

        VertexPosMap        M_original_pos_map;
        std::mt19937_64     randomizer;
        float               last_x1 = 0, last_x2 = 0, last_x3 = 0;
//        static std::mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31,
//        0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253>
    };
    
}

#endif /* defined(__MeshEditE__test_random_rotation__) */
