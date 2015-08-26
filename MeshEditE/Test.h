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
#define ROUNDER       100000.0
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
CGLA::Vec3d     mul_3D_dir                      ( const CGLA::Mat4x4d &t, const CGLA::Vec3d &dir );

double          sphere_intersects               ( const CGLA::Vec3d &c1, double r1,
                                                  const CGLA::Vec3d &c2, double r2 );

inline double truncateDouble( double d ){
   return ( std::floor( d * ROUNDER ) / ROUNDER );
}

inline double truncateDouble( double d, size_t digits ){
    size_t rounder = 10;
    for( size_t i = 1; i < digits; ++i ){ rounder *= 10;}
    double rounder_f = static_cast<double>( rounder );
    return ( std::floor( d * rounder_f ) / rounder_f );
}

inline double truncateDouble2( double d ){ return truncateDouble( d, 2 ); }
inline double truncateDouble3( double d ){ return truncateDouble( d, 3 ); }
inline double truncateDouble4( double d ){ return truncateDouble( d, 4 ); }
inline double truncateDouble5( double d ){ return truncateDouble( d, 5 ); }
inline double truncateDouble6( double d ){ return truncateDouble( d, 6 ); }
inline double truncateDouble7( double d ){ return truncateDouble( d, 7 ); }
inline void truncateVec3d( CGLA::Vec3d& v ){
    v[0] = truncateDouble4( v[0] );
    v[1] = truncateDouble4( v[1] );
    v[2] = truncateDouble4( v[2] );
}

inline void truncateMat4x4d( CGLA::Mat4x4d& M ){
    M[0][0] = truncateDouble4( M[0][0] );
    M[0][1] = truncateDouble4( M[0][1] );
    M[0][2] = truncateDouble4( M[0][2] );
    M[0][3] = truncateDouble4( M[0][3] );

    M[1][0] = truncateDouble4( M[1][0] );
    M[1][1] = truncateDouble4( M[1][1] );
    M[1][2] = truncateDouble4( M[1][2] );
    M[1][3] = truncateDouble4( M[1][3] );
    
    M[2][0] = truncateDouble4( M[2][0] );
    M[2][1] = truncateDouble4( M[2][1] );
    M[2][2] = truncateDouble4( M[2][2] );
    M[2][3] = truncateDouble4( M[2][3] );

    M[3][0] = truncateDouble4( M[3][0] );
    M[3][1] = truncateDouble4( M[3][1] );
    M[3][2] = truncateDouble4( M[3][2] );
    M[3][3] = truncateDouble4( M[3][3] );
}


inline void checkVec3( CGLA::Vec3d v ){
    assert( isfinite( v[0] )); assert( isfinite( v[1] )); assert( isfinite( v[2] ));
}

inline void checkMat4( CGLA::Mat4x4d M ){
    assert( isfinite( M[0][0] )); assert( isfinite( M[0][1] )); assert( isfinite( M[0][2] ));     assert( isfinite( M[0][3] ));
    assert( isfinite( M[1][0] )); assert( isfinite( M[1][1] )); assert( isfinite( M[1][2] ));     assert( isfinite( M[1][3] ));
    assert( isfinite( M[2][0] )); assert( isfinite( M[2][1] )); assert( isfinite( M[2][2] ));     assert( isfinite( M[2][3] ));
    assert( isfinite( M[3][0] )); assert( isfinite( M[3][1] )); assert( isfinite( M[3][2] ));     assert( isfinite( M[3][3] ));
}

inline void mat4Copy( const CGLA::Mat4x4d& source, CGLA::Mat4x4d& dest ){
 
    checkMat4( source );
    dest[0][0] = source[0][0];  dest[0][1] = source[0][1]; dest[0][2] = source[0][2];  dest[0][3] = source[0][3];
    dest[1][0] = source[1][0];  dest[1][1] = source[1][1]; dest[1][2] = source[1][2];  dest[1][3] = source[1][3];
    dest[2][0] = source[2][0];  dest[2][1] = source[2][1]; dest[2][2] = source[2][2];  dest[2][3] = source[2][3];
    dest[3][0] = source[3][0];  dest[3][1] = source[3][1]; dest[3][2] = source[3][2];  dest[3][3] = source[3][3];
    checkMat4( dest );
}

inline template<typename T>
T            in_range                        ( T value, T low, T high )
{
    if( value < low  ) return low;
    if( value > high ) return high;
    return value;
}




#endif /* defined(__MeshEditE__Test__) */
