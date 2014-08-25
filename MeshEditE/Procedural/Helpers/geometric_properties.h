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
#include "polarize.h"


namespace Procedural{
    namespace Geometry{

typedef std::pair< int, double > DistanceMetrics;

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
double              angle                                   ( CGLA::Vec3d l, CGLA::Vec3d r );
// must have calld label junction on edge_info
void                distance_from_poles                     ( HMesh::Manifold& m,
                                                              HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info,
                                                              HMesh::VertexAttributeVector<DistanceMetrics> &distances );
void                distance_from_junctions                 ( HMesh::Manifold& m,
                                                              HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info,
                                                              HMesh::VertexAttributeVector<DistanceMetrics> &distances );
void                distance_from_poles_and_junctions       ( HMesh::Manifold& m,
                                                              HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info,
                                                              HMesh::VertexAttributeVector<DistanceMetrics> &distances );




}}


#endif /* defined(__MeshEditE__geometric_properties__) */
