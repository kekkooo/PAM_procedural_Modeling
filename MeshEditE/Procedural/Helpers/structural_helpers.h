//
//  structural_helpers.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__structural_helpers__
#define __MeshEditE__structural_helpers__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>
#include "polarize.h"


namespace Procedural{
    namespace Structure{


HMesh::HalfEdgeID   find_half_edge      ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1);

HMesh::HalfEdgeID   find_rib_edge       ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1,
                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info);
HMesh::HalfEdgeID   find_spine_edge     ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1,
                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info);
}}

#endif /* defined(__MeshEditE__structural_helpers__) */
