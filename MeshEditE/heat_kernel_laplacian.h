//
//  heat_kernel_laplacian.h
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 05/02/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__heat_kernel_laplacian__
#define __MeshEditE__heat_kernel_laplacian__

#include <iostream>
#include <GEL/HMesh/Manifold.h>

HMesh::VertexAttributeVector<double> heat_kernel_laplacian(HMesh::Manifold& m, double R = 2.0);


#endif /* defined(__MeshEditE__heat_kernel_laplacian__) */
