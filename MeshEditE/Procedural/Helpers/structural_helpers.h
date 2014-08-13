//
//  structural_helpers.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__structural_helpers__
#define __MeshEditE__structural_helpers__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>
#include "polarize.h"



namespace Procedural{
    namespace Structure{

enum EndType{ Invalid, Junction, Pole };
struct RegionBoundaries
{
    HMesh::HalfEdgeID   FirstSpineEdge;
    EndType             FirstEndType;
    HMesh::HalfEdgeID   SecondSpineEdge;
    EndType             SecondEndType;
    
    RegionBoundaries( HMesh::HalfEdgeID FirstEdge = HMesh::InvalidHalfEdgeID,
                      EndType FirstEnd = EndType::Invalid,
                      HMesh::HalfEdgeID SecondEdge = HMesh::InvalidHalfEdgeID,
                     EndType SecondEnd = EndType::Invalid ):
    FirstSpineEdge(FirstEdge),   FirstEndType(FirstEnd),
    SecondSpineEdge(SecondEdge), SecondEndType(SecondEnd) { }
                     
};
        

HMesh::HalfEdgeID   find_half_edge      ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1);

HMesh::HalfEdgeID   find_rib_edge       ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1,
                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info);
HMesh::HalfEdgeID   find_spine_edge     ( HMesh::Manifold& m, HMesh::VertexID v0, HMesh::VertexID v1,
                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info);
        
RegionBoundaries    FollowSpines        ( HMesh::Manifold& m, HMesh::VertexID,
                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info );
// ALERT : IS NOT THAT EASY. YOU NEED TO BUILD THE GRAPH IN ORDER TO FIND LOOPS
HMesh::HalfEdgeID   BranchCutDirection  ( HMesh::Manifold& m, HMesh::VertexID,
                                          HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info );
}}

#endif /* defined(__MeshEditE__structural_helpers__) */
