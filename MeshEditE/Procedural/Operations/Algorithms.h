//
//  Algorithms.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Algorithms__
#define __MeshEditE__Algorithms__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <polarize.h>

namespace Procedural{
    namespace Operations{
        namespace Algorithms{

    void selected_vertices_cotangent_weights_laplacian  ( HMesh::Manifold&m, std::vector< HMesh::VertexID > vs );
    void selected_vertices_inverse_distance_laplacian   ( HMesh::Manifold&m, std::vector< HMesh::VertexID > vs );
    void inverse_distance_laplacian_smoothing           ( HMesh::Manifold& am );
    void cotangent_weights_laplacian_smoothing          ( HMesh::Manifold& am );
    void along_spines                                   ( HMesh::Manifold& m,
                                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> &edge_info );
//    void along_ribs                                     ( HMesh::Manifold& m,
//                                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info );

}}}

#endif /* defined(__MeshEditE__Algorithms__) */
