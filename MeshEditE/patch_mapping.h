//
//  patch_mapping.h
//  MeshEditE
//
//  Created by Francesco Usai on 06/10/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__patch_mapping__
#define __MeshEditE__patch_mapping__

#include <stdio.h>
#include <GEL/HMesh/Manifold.h>

void build_patches( HMesh::Manifold& m, HMesh::FaceAttributeVector<int> &face_to_patch );
bool is_corner ( HMesh::Manifold& m, HMesh::VertexID v, HMesh::HalfEdgeAttributeVector<int> &layout );
inline bool is_corner ( HMesh::VertexID v, std::vector<HMesh::VertexID> &corners )
{
    return ( std::find( corners.begin(), corners.end(), v ) != corners.end() );
}
#endif /* defined(__MeshEditE__patch_mapping__) */
