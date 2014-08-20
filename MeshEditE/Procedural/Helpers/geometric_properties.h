//
//  geometric_properties.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__geometric_properties__
#define __MeshEditE__geometric_properties__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>

namespace Procedural{
    namespace Geometry{


CGLA::Vec3d         face_normal                             ( HMesh::Manifold& m, HMesh::FaceID f       );
CGLA::Vec3d         vertex_normal                           ( HMesh::Manifold& m, HMesh::VertexID v     );
CGLA::Vec3d         ring_barycenter                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h   );
CGLA::Vec3d         ring_barycenter                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                             std::vector< HMesh::VertexID > &vertices   );

CGLA::Vec3d         ring_barycenter                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                             std::vector< HMesh::VertexID > &vertices,
                                                             std::vector< HMesh::HalfEdgeID > &halfedges);

void                ring_vertices_and_halfedges             ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                             std::vector< HMesh::VertexID > &vertices,
                                                             std::vector< HMesh::HalfEdgeID > &halfedges);

int                 valence                                 ( HMesh::Manifold& m, HMesh::VertexID v     );
bool                is_singularity                          ( HMesh::Manifold& m, HMesh::VertexID v     );
bool                is_2_neighbor_of_pole                   ( HMesh::Manifold& m, HMesh::VertexID v     );

}}


#endif /* defined(__MeshEditE__geometric_properties__) */
