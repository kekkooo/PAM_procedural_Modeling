//
//  basic_shapes.h
//  MeshEditE
//
//  Created by Francesco Usai on 06/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__basic_shapes__
#define __MeshEditE__basic_shapes__

#include <iostream>
#include <GEL/HMesh/Manifold.h>

namespace Procedural
{
namespace Operations
{
    void    create_triangle_cube    ( HMesh::Manifold& m, double edge_length    );
    void    create_quads_cube       ( HMesh::Manifold& m, double edge_length    );
    void    create_basic_PAM        ( HMesh::Manifold& m, double edge_length    );
    void    create_PAM_box          ( HMesh::Manifold& m, double l1, double l2, double l3 );
}
}

#endif /* defined(__MeshEditE__basic_shapes__) */
