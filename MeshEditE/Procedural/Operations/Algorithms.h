//
//  Algorithms.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Algorithms__
#define __MeshEditE__Algorithms__

#include <iostream>
#include <GEL/HMesh/Manifold.h>

namespace Procedural{
    namespace Operations{
        namespace Algorithms{

    void        inverse_distance_laplacian_smoothing    ( HMesh::Manifold& am );
    void        cotangent_weights_laplacian_smoothing   ( HMesh::Manifold& am );

}}}

#endif /* defined(__MeshEditE__Algorithms__) */
