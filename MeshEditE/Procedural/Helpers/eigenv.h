//
//  eigenv.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/07/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef MeshEditE_eigenv_h
#define MeshEditE_eigenv_h

#include <vector>
#include <GEL/CGLA/Vec3d.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace Procedural{
    namespace Geometry{
        
        struct PrincipalDirection{
            CGLA::Vec3d dir;
            double      magnitude;
        };
        
        struct PCAResult{
            PrincipalDirection major;
            PrincipalDirection minor;
            PrincipalDirection other;
        };
        
        static Eigen::Matrix3d outerProduct( const CGLA::Vec3d& p0, const CGLA::Vec3d& p1 ){
            CGLA::Vec3d r1, r2, r3;
            Eigen::Matrix3d a = Eigen::Matrix3d::Zero();
            r1 = p1*p0[0];
            r2 = p1*p0[1];
            r3 = p1*p0[2];
            
            a << r1[0] , r1[1] , r1[2],
            r2[0] , r2[1] , r2[2],
            r3[0] , r3[1] , r3[2];
            
            //        std::cout << p0 << std::endl << p1 << std::endl << a << std::endl;
            
            return a;
        }
    }
    
    static Geometry::PCAResult EigenVectors( const std::vector<CGLA::Vec3d>& points ){
        
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        
        Vec3d mean(0);
        for( const CGLA::Vec3d& p : points ){
            mean[0] += p[0];
            mean[1] += p[1];
            mean[2] += p[2];
        }
        double no_points = static_cast<double>( points.size( ));
        mean /= no_points;
        for( const CGLA::Vec3d& po : points ){
            Vec3d p = po - mean;
            Eigen::Matrix3d a = Geometry::outerProduct( p, p );
            cov = cov + a;
        }
        
        cov = cov / ( no_points - 1.0 );
        
        std::cout << " Covariance Matrix " << std::endl << cov << std::endl;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(cov);
        
        if( eigensolver.info() != Eigen::Success ){
            std::cout << "cannot decompose covariance matrix "<< std::endl;
            assert(false);
        }
        else{
            std::cout << "eigenvectors are : " << eigensolver.eigenvectors() << std::endl
            << "eigenvalues  are : " << eigensolver.eigenvalues() << std::endl;
        }
        
        Geometry::PrincipalDirection major, minor, other;
        Geometry::PCAResult p;
        other.dir[0] = eigensolver.eigenvectors()( 0, 0 );
        other.dir[1] = eigensolver.eigenvectors()( 1, 0 );
        other.dir[2] = eigensolver.eigenvectors()( 2, 0 );
        other.magnitude = eigensolver.eigenvalues()( 0) ;
        
        minor.dir[0] = eigensolver.eigenvectors()( 0, 1 );
        minor.dir[1] = eigensolver.eigenvectors()( 1, 1 );
        minor.dir[2] = eigensolver.eigenvectors()( 2, 1 );
        minor.magnitude = eigensolver.eigenvalues()( 1 );

        
        major.dir[0] = eigensolver.eigenvectors()( 0, 2 );
        major.dir[1] = eigensolver.eigenvectors()( 1, 2 );
        major.dir[2] = eigensolver.eigenvectors()( 2, 2 );
        major.magnitude = eigensolver.eigenvalues()( 2 );
        
        std::cout << " other " << std::endl << other.dir << ", " << other.magnitude << std::endl;
        std::cout << " minor " << std::endl << minor.dir << ", " << minor.magnitude << std::endl;
        std::cout << " major " << std::endl << major.dir << ", " << major.magnitude << std::endl;
        
        p.other = other;
        p.minor = minor;
        p.major = major;
        
        return p;
    }
    
}

#endif
