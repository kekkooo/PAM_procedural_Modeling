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


// Geometric Structure Utilities : structural_helpers.h



// other stuff                      : other.h
bool                test_ring_barycenter                    ( HMesh::Manifold& m, HMesh::HalfEdgeID h );
CGLA::Vec3d         simple_random_direction                 ( HMesh::Manifold& m, HMesh::VertexID   v, double face_normal_weight );
CGLA::Vec3d         alt_simple_random_direction             ( HMesh::Manifold& m, HMesh::VertexID   v );
CGLA::Vec3f         color_ramp                              ( int value, int max );
CGLA::Vec3f         color_ramp2                             ( int value, int max );
void                linspace                                ( double min, double max, int num,
                                                             std::vector<double> &values);


#endif /* defined(__MeshEditE__Test__) */
