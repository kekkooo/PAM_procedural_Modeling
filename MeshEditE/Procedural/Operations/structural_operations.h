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
#include <polarize.h>



namespace Procedural{
    namespace Operations{
        namespace Structural{

HMesh::VertexID add_branch     ( HMesh::Manifold& m, HMesh::VertexID, int size, HMesh::VertexAttributeVector<int> &ring );
void            cut_branch     ( HMesh::Manifold& m, HMesh::HalfEdgeID h );
void            cut_branch     ( HMesh::Manifold& m, HMesh::HalfEdgeID h, HMesh::VertexID pole);
void            cut_branch     ( HMesh::Manifold& m, HMesh::VertexID v, HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info );
void            remove_branch  ( HMesh::Manifold& m, HMesh::VertexID pole, HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info );
void            glue_poles_with_valence_equalization
                               ( HMesh::Manifold& m, HMesh::VertexID pole1, HMesh::VertexID pole2 );
bool            glue_poles     ( HMesh::Manifold &m, HMesh::VertexID pole1, HMesh::VertexID pole2 );

}}}
#endif /* defined(__MeshEditE__structural_opeations__) */
