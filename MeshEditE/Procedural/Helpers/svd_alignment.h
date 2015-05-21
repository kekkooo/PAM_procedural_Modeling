//
//  svd_alignment.h
//  MeshEditE
//
//  Created by Francesco Usai on 23/09/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__svd_alignment__
#define __MeshEditE__svd_alignment__

#include <iostream>
#include <GEL/LinAlg/Matrix.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Mat4x4d.h>

namespace Procedural {
    namespace Geometry {

        double Det( LinAlg::CMatrix m );
        
        void svd_rigid_motion( std::vector<CGLA::Vec3d> &P, std::vector<CGLA::Vec3d> &Q,
                               CGLA::Mat4x4d &rot, CGLA::Mat4x4d &translation );

        /// calculate the rigid transformation that aligns P w.r.t Q
        void svd_rigid_motion( HMesh::Manifold& m1, std::vector< HMesh::VertexID > &P,
                               HMesh::Manifold& m2, std::vector< HMesh::VertexID > &Q,
                               CGLA::Mat4x4d &rot, CGLA::Mat4x4d &translation );
}}

#endif /* defined(__MeshEditE__svd_alignment__) */
