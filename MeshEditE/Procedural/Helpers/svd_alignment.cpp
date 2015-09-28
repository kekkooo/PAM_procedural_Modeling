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

#include <Eigen/Dense>
#include <Eigen/SVD>


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
        
        
        void svd_rigid_motion( const vector<Vec3d> &P, const vector<Vec3d> &Q, const vector<double> &weights, Mat4x4d &rot, Mat4x4d &translation ){
            assert( P.size() == Q.size( ));
            assert( P.size() == weights.size( ));
            int n = P.size();
            // assume weights are always 1
            // 1) find the centroid of the two point sets
            Vec3d p( 0.0, 0.0, 0.0 ),
                  q( 0.0, 0.0, 0.0 );
            
            for( Vec3d Pi : P ) { p += Pi; }
            for( Vec3d Qi : Q ) { q += Qi; }
            
            p /= n;
            q /= n;
            
            // 2) compute the centered vectors
            vector< Vec3d > xs, ys;
            for( const Vec3d& Pi : P ) { xs.push_back( Pi - p ); }
            for( const Vec3d& Qi : Q ) { ys.push_back( Qi - q ); }
            
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
                
                W.set( i, i, weights[i] );
            }
            CMatrix S = X * W * Y_t;
            // it should be 3 * 3
            assert( S.Rows() == 3 && S.Cols() == 3 );
            cout << S.Rows() << " # " << S.Cols() << endl;
            // ask if V is transposed or not, but I guess it is
            // given the signature of the method it seems to be not.
            CMatrix U, Sigma, V;
            SVD( S, U, Sigma, V);
            
            CMatrix U_t = U.Transposed(),
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
            
            Vec3d   t   = q - rot.mul_3D_point( p );
            translation = CGLA::translation_Mat4x4d( t );
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
            for( VertexID& vid : P ) { xs.push_back( m1.pos(vid) - p ); }
            for( VertexID& vid : Q ) { ys.push_back( m2.pos(vid) - q ); }
            
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
        
        
        void svd_6d_rigid_motion( const std::vector<CGLA::Vec3d> &P, const std::vector<CGLA::Vec3d> &Q,
                                  const std::vector<CGLA::Vec3d> &Pn, const std::vector<CGLA::Vec3d> &Qn,
                                  CGLA::Mat4x4d &rot, CGLA::Mat4x4d &translation ){
            assert( P.size() == Q.size( ));
            assert( P.size() == Pn.size( ));
            assert( P.size() == Qn.size( ));
            size_t n = P.size();
            
            // convert point and normal to a 6D vector
            // and calculate the centroid of the 6D points
            typedef Eigen::VectorXd vec6D;
            std::vector<Eigen::VectorXd> ps, qs;
            Eigen::VectorXd cps(6), cqs(6);
            cps.setZero();
            cqs.setZero();
//            cps << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
//            cqs << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
            
            for( int i = 0; i < P.size(); ++i ){
                Eigen::VectorXd vp(6), vq(6);
                const CGLA::Vec3d& pi   = P.at( i );
                const CGLA::Vec3d& pni  = Pn.at( i );
                const CGLA::Vec3d& qi   = Q.at( i );
                const CGLA::Vec3d& qni  = Qn.at( i );
                
                // the normal of vq should be opposed, since we want them to be aligned and opposite
                vp << pi[0], pi[1], pi[2], pni[0], pni[1], pni[2];
//                vq << qi[0], qi[1], qi[2], qni[0], qni[1], qni[2];
                vq << qi[0], qi[1], qi[2], -qni[0], -qni[1], -qni[2];
                
                ps.push_back( vp );
                qs.push_back( vq );
                
                cps = cps + vp;
                
                cqs = cqs + vq;

            }
            
            cps /= static_cast<double>( n );
            cqs /= static_cast<double>( n );
            
//            Vec3d cp_normal( cps( 3 ), cps( 4 ), cps( 5 ));
//            cp_normal.normalize();
//            cps( 3 ) = cp_normal[0];
//            cps( 4 ) = cp_normal[1];
//            cps( 5 ) = cp_normal[2];
//
//            Vec3d cq_normal( cqs( 3 ), cqs( 4 ), cqs( 5 ));
//            cq_normal.normalize();
//            cqs( 3 ) = cq_normal[0];
//            cqs( 4 ) = cq_normal[1];
//            cqs( 5 ) = cq_normal[2];

            
            // 2) compute the centered vectors
            std::vector<Eigen::VectorXd> xs, ys;
            for( const auto& pi : ps) {
                Eigen::VectorXd xi = pi - cps;

//                Vec3d xi_normal( xi( 3 ), xi( 4 ), xi( 5 ));
//                xi_normal.normalize();
//                xi( 3 ) = xi_normal[0];
//                xi( 4 ) = xi_normal[1];
//                xi( 5 ) = xi_normal[2];

                xs.push_back( xi );
            }
            for( const auto& qi : qs) {
                Eigen::VectorXd yi = qi - cqs;
                
//                Vec3d yi_normal( yi( 3 ), yi( 4 ), yi( 5 ));
//                yi_normal.normalize();
//                yi( 3 ) = yi_normal[0];
//                yi( 4 ) = yi_normal[1];
//                yi( 5 ) = yi_normal[2];

                ys.push_back( yi );
            }
            
            // 3) compute D * D covariance matrix S = XWY_t
            // X and Y are 6 * n matrices that have x_i and y_i as their columns
            // and W is diag( w0, w1, ..., wn ) ( in my case is I )

            typedef Eigen::Matrix<double, 6,6> Mat6x6d;
            Mat6x6d X, Y_t;
            Mat6x6d W = Mat6x6d::Identity();
            // set the X, Y matrices
            for( int i = 0; i < n; ++i )
            {
                X( 0, i ) = xs[i]( 0 );
                X( 1, i ) = xs[i]( 1 );
                X( 2, i ) = xs[i]( 2 );
                X( 3, i ) = xs[i]( 3 );
                X( 4, i ) = xs[i]( 4 );
                X( 5, i ) = xs[i]( 5 );
                
                Y_t( i, 0 ) = ys[i]( 0 );
                Y_t( i, 1 ) = ys[i]( 1 );
                Y_t( i, 2 ) = ys[i]( 2 );
                Y_t( i, 3 ) = ys[i]( 3 );
                Y_t( i, 4 ) = ys[i]( 4 );
                Y_t( i, 5 ) = ys[i]( 5 );
            }

            
            Mat6x6d S = X * W * Y_t;
            // it should be 6 * 6

            Mat6x6d U, Sigma, V;
            Mat6x6d U_t, V_t;
            
            Eigen::JacobiSVD<Mat6x6d> svd( S, Eigen::ComputeFullU | Eigen::ComputeFullV );
            
            U = svd.matrixU();
            V = svd.matrixV();
            U_t = U.transpose();
            V_t = V.transpose();

            Mat6x6d V_U_t = V * U_t;
            
            double det_V_U_t = V_U_t.determinant();
            
            Mat6x6d M = Mat6x6d::Identity();
            M( 5, 5 ) = det_V_U_t;
            Mat6x6d R = V * M * U_t;
            
            
            // go back to CGLA
            rot = CGLA::Mat4x4d( 0.0 );
            for( int i = 0; i < 3; ++i )
            {
                for( int j = 0; j < 3; ++j)
                {
                    rot[i][j] = R( i, j );
                }
            }
            rot[3][3] = 1;
            
            Vec3d CGLA_cqs( cqs( 0 ), cqs( 1 ), cqs( 2 ));
            Vec3d CGLA_cps( cps( 0 ), cps( 1 ), cps( 2 ));
            
            Vec3d   t   = CGLA_cqs - rot.mul_3D_point( CGLA_cps );
            translation = CGLA::translation_Mat4x4d( t );
        }


}}