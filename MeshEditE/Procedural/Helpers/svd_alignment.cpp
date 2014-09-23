//
//  svd_alignment.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 23/09/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "svd_alignment.h"
#include <GEL/CGLA/Vec3d.h>
#include <GEL/LinAlg/Vector.h>
#include <GEL/LinAlg/LapackFunc.h>

using namespace HMesh;
using namespace LinAlg;
using namespace std;
using namespace CGLA;

namespace Procedural {
    namespace Geometry {
        
        double Det( CMatrix m )
        {
            // the matrix should be squared
            assert(m.Rows() == m.Cols());
            
            double det =
            (
             m.get(0, 0) * m.get(1, 1) * m.get(2, 2)
             + m.get(0, 1) * m.get(1, 2) * m.get(2, 0)
             + m.get(0, 2) * m.get(1, 0) * m.get(2, 1)
             )
            -
            (
             m.get(0, 2) * m.get(1, 1) * m.get(2, 0)
             + m.get(0, 1) * m.get(1, 0) * m.get(2, 2)
             + m.get(0, 0) * m.get(1, 2) * m.get(2, 1)
             );
            return det;
        }
        
        /// calculate the rigid transformation that aligns P w.r.t Q
        void svd_rigid_motion( Manifold& m1, vector< VertexID > &P,
                               Manifold& m2, vector< VertexID > &Q,
                               Mat4x4d &rot, Mat4x4d &translation )
        {
            assert( P.size() == Q.size( ));
            int n = P.size();
            // assume weights are always 1
            // 1) find the centroid of the two point sets
            Vec3d p( 0.0, 0.0, 0.0 ),
                  q( 0.0, 0.0, 0.0 );
            
            for( VertexID vid : P ) { p += m1.pos( vid ); }
            for( VertexID vid : Q ) { q += m2.pos( vid ); }
            
            p /= n;
            q /= n;
            
            // 2) compute the centered vectors
            vector< Vec3d > xs, ys;
            for( VertexID vid : P ) { xs.push_back( m1.pos(vid) - p ); }
            for( VertexID vid : Q ) { ys.push_back( m2.pos(vid) - q ); }
            
            // 3) compute D * D covariance matrix S = XWY_t
            // X and Y are 3 * n matrices that have x_i and y_i as their columns
            // and W is diag( w0, w1, ..., wn ) ( in my case is I )
            LinAlg::CMatrix X( 3, n, 0.0 ), W( n, n, 0.0 ), Y_t( n, 3, 0.0 );
            for( int i = 0; i < n; ++i )
            {
                X.set( 0, i, xs[i][0] );
                X.set( 1, i, xs[i][1] );
                X.set( 2, i, xs[i][2] );
                
                Y_t.set( i, 0, ys[i][0] );
                Y_t.set( i, 1, ys[i][1] );
                Y_t.set( i, 2, ys[i][2] );
                
                W.set( i, i, 1.0 );
            }
            CMatrix S = X * W * Y_t;
            // it should be 3 * 3
            assert( S.Rows() == 3 && S.Cols() == 3 );
            cout << S.Rows() << " # " << S.Cols() << endl;
            // ask if V is transposed or not, but I guess it is
            // given the signature of the method it seems to be not.
            CMatrix U, Sigma, V;
            SVD( S, U, Sigma, V);
            
            CMatrix U_t         = U.Transposed(),
                    V_t         = V.Transposed(),
                    V_U_t       = ( V * U_t );
            double  det_V_U_t   = Det( V_U_t );
            CMatrix M( 3, 3, 0 );
            M.set( 0, 0, 1.0 );
            M.set( 1, 1, 1.0 );
            M.set( 2, 2, det_V_U_t );
            
            CMatrix R = V * M * U_t;
            assert( R.Rows() == 3 && R.Cols() == 3 );
            // go back to CGLA
            rot = CGLA::Mat4x4d( 0.0 );
            for( int i = 0; i < 3; ++i )
            {
                for( int j = 0; j < 3; ++j)
                {
                    rot[i][j] = R.get( i, j );
                }
            }
            rot[3][3] = 1;
            
            Vec3d   t           = q - rot.mul_3D_point( p );
                    translation = CGLA::translation_Mat4x4d( t );
        }

}}