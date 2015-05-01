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
#include <GEL/CGLA/Mat4x4d.h>
#include <GEL/CGLA/Mat4x4f.h>


// Geometric Structure Utilities : structural_helpers.h
#define GEO_EPS     0.00000001
#define ARITH_EPS   0.000001


// other stuff                      : other.h
bool            test_ring_barycenter            ( HMesh::Manifold& m, HMesh::HalfEdgeID h );
CGLA::Vec3d     simple_random_direction         ( HMesh::Manifold& m, HMesh::VertexID   v, double face_normal_weight );
CGLA::Vec3d     alt_simple_random_direction     ( HMesh::Manifold& m, HMesh::VertexID   v );
CGLA::Vec3f     color_ramp                      ( int value, int max );
CGLA::Vec3f     color_ramp2                     ( int value, int max );
void            linspace                        ( double min, double max, int num,
                                                      std::vector<double> &values);
CGLA::Mat4x4d   get_rotation_mat4d              ( CGLA::Vec3d axis, double cosine );
CGLA::Mat4x4d   get_alignment_for_2_vectors     ( CGLA::Vec3d v1, CGLA::Vec3d v2 );
void            alt_glue_poles                  ( HMesh::Manifold& mani, HMesh::VertexID vid0, HMesh::VertexID vid1);
void            bezier                          ( CGLA::Vec3d p0, CGLA::Vec3d p1, CGLA::Vec3d p2, int n, std::vector< CGLA::Vec3d >& points );
void            save_colored_obj                ( HMesh::Manifold &m, std::string &path );
void            save_intermediate_result        ( HMesh::Manifold &m, const std::string &path, int step_number );
// takes two manifolds and two poles and returns the number of ccw nexts should be done
// on p2, in order to minimize the angle distance between the two triangle fans
int             get_starter_offset              ( HMesh::Manifold &m1, HMesh::VertexID p1,
                                                  HMesh::Manifold &m2, HMesh::VertexID p2 );
void            bridge_pole_one_rings           ( HMesh::Manifold& mani, HMesh::VertexID vid0, HMesh::VertexID vid1);
void            buildReflectionMatrix           ( CGLA::Vec3d& planeNormal, CGLA::Mat4x4d &T );
CGLA::Mat4x4d   alt_get_alignment_for_2_vectors ( CGLA::Vec3d v1, CGLA::Vec3d v2 );
CGLA::Mat4x4d   Mat4x4d_to_float                ( CGLA::Mat4x4f &m );

inline template<typename T>
T            in_range                        ( T value, T low, T high )
{
    if( value < low  ) return low;
    if( value > high ) return high;
    return value;
}



#endif /* defined(__MeshEditE__Test__) */
