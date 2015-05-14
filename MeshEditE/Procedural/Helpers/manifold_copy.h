//
//  manifold_copy.h
//  MeshEditE
//
//  Created by Francesco Usai on 06/11/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__manifold_copy__
#define __MeshEditE__manifold_copy__

#include <stdio.h>
#include <GEL/HMesh/Manifold.h>
#include <set>

namespace Procedural{
    namespace Helpers{
        
        void add_manifold( HMesh::Manifold &m, HMesh::Manifold &other,
                           std::set<HMesh::VertexID> &m_IDs_in_result,
                           std::set<HMesh::VertexID> & other_IDs_in_result );
        
        void add_manifold( HMesh::Manifold &m, HMesh::Manifold &other,
                           std::set<HMesh::VertexID> &m_IDs_in_result,
                           std::set<HMesh::VertexID> & other_IDs_in_result,
                           HMesh::IDRemap &remap );

        
        void test_delete ( HMesh::Manifold &m, std::set<HMesh::VertexID> &to_delete );
        
}}

#endif /* defined(__MeshEditE__manifold_copy__) */
