//
//  Test.h
//  MeshEditE
//
//  Created by Francesco Usai on 25/07/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Test__
#define __MeshEditE__Test__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <tuple>
#include "polarize.h"

// DATA TYPES
typedef std::vector< std::tuple< HMesh::VertexID, double, CGLA::Vec3d > >  DistanceDirectionVector;
enum VertexType
{
    REGULAR     = 0x01,
    SINGULAR    = 0x02,
    POLE        = 0x04
};
// JUST BASIC STUFF -> into test.
void                create_triangle_cube                    ( HMesh::Manifold& m, double edge_length    );
void                create_quads_cube                       ( HMesh::Manifold& m, double edge_length    );
void                create_basic_PAM                        ( HMesh::Manifold& m, double edge_length    );

// SMOOTHING ALGORITHMS -> test.h
void                inverse_distance_laplacian_smoothing    ( HMesh::Manifold& am                       );
void                cotangent_weights_laplacian_smoothing   ( HMesh::Manifold& am                       );

// OPERATIONS on GEOMETRY : -> geometric_operations.h
// h is a spine edge
std::vector< HMesh::VertexID >
                    split_ring_of_quads                     ( HMesh::Manifold& m, HMesh::HalfEdgeID h   );
// h is a rib edge
std::vector< HMesh::VertexID >
                    split_from_pole_to_pole                 ( HMesh::Manifold& m, HMesh::HalfEdgeID h   );


void                scale_ring_radius                       ( HMesh::Manifold& m, HMesh::HalfEdgeID h, double ratio );
void                set_ring_radius                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h, double radius );
void                set_ring_radius                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h, std::vector<double> radii );
void                set_ring_radius                         ( HMesh::Manifold& m, CGLA::Vec3d barycenter,
                                                              DistanceDirectionVector ddv );

void                flatten_pole                            ( HMesh::Manifold& m, HMesh::VertexID pole );
void                extrude_pole                            ( HMesh::Manifold& m, HMesh::VertexID v,
                                                             double length = -1.0, bool pole_add_quads = false,
                                                             double radius_preserving_ratio = 1.0      );
void                extrude_pole                            ( HMesh::Manifold&m, HMesh::VertexID v, CGLA::Vec3d direction,
                                                              bool add_quads, double scaling );

void                move_vertex                             ( HMesh::Manifold& m, HMesh::VertexID v,
                                                             double length = -1.0, bool pole_add_quads = false,
                                                             double radius_preserving_ratio = 1.0      );

void                add_noise                               ( HMesh::Manifold& m, int vertex_type_flags, double ratiom, double cutoff );

// OPERATIONS on the STRUCTURE  -> structural_operations.h
void                cut_branch                              ( HMesh::Manifold& m, HMesh::HalfEdgeID h );
void                cut_branch                              ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                              HMesh::VertexID pole                      );
void                add_ring_around_pole                    ( HMesh::Manifold& m, HMesh::VertexID pole,
                                                             double scaling );



// Geometric Properties Utilities : geometric_properties.h
CGLA::Vec3d         face_normal                             ( HMesh::Manifold& m, HMesh::FaceID f       );
CGLA::Vec3d         vertex_normal                           ( HMesh::Manifold& m, HMesh::VertexID v     );
CGLA::Vec3d         ring_barycenter                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h   );
CGLA::Vec3d         ring_barycenter                         ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                              std::vector< HMesh::VertexID > &vertices  );
int                 valence                                 ( HMesh::Manifold& m, HMesh::VertexID v     );
bool                is_singularity                          ( HMesh::Manifold& m, HMesh::VertexID v     );

// Geometric Structure Utilities : structural_helpers.h
HMesh::HalfEdgeID   find_half_edge                          ( HMesh::Manifold& m, HMesh::VertexID v0,
                                                              HMesh::VertexID v1                        );

HMesh::HalfEdgeID   find_rib_edge                           ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1,
                                                              HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info);
HMesh::HalfEdgeID   find_spine_edge                         ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1,
                                                              HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info);


// other stuff                      : other.h
bool                test_ring_barycenter                    ( HMesh::Manifold& m, HMesh::HalfEdgeID h );
CGLA::Vec3d         simple_random_direction                 ( HMesh::Manifold& m, HMesh::VertexID   v, double face_normal_weight );
CGLA::Vec3d         alt_simple_random_direction             ( HMesh::Manifold& m, HMesh::VertexID   v );


#endif /* defined(__MeshEditE__Test__) */
