//
//  structural_opeations.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__structural_opeations__
#define __MeshEditE__structural_opeations__

#include <iostream>
#include <GEL/HMesh/Manifold.h>



namespace Procedural{
    namespace Operations{
        namespace Structural{


void                cut_branch                              ( HMesh::Manifold& m, HMesh::HalfEdgeID h );
void                cut_branch                              ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                              HMesh::VertexID pole                      );
            // this should be placed into geometric operationss

}}}
#endif /* defined(__MeshEditE__structural_opeations__) */
