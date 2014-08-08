//
//  structural_helpers.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "structural_helpers.h"

using namespace HMesh;

namespace Procedural{
    namespace Structure{

HalfEdgeID find_half_edge( Manifold& m, VertexID v0, VertexID v1 )
{
    Walker w = m.walker(v0);
    HalfEdgeID found = InvalidHalfEdgeID;
    for( ; !w.full_circle() && found == InvalidHalfEdgeID; w = w.circulate_vertex_ccw())
    {
        if( w.vertex() == v1 )
        {
            found = w.halfedge();
        }
    }
    return found;
}

HalfEdgeID find_rib_edge( Manifold& m, VertexID v0, VertexID v1,
                         HalfEdgeAttributeVector<EdgeInfo> edge_info)
{
    HalfEdgeID found = find_half_edge( m, v0, v1 );
    if ( found == InvalidHalfEdgeID || !edge_info[found].is_rib( ))
    {
        return InvalidHalfEdgeID;
    }
    else
    {
        return found;
    }
}

HalfEdgeID find_spine_edge( Manifold& m, VertexID v0, VertexID v1,
                           HalfEdgeAttributeVector<EdgeInfo> edge_info)
{
    HalfEdgeID found = find_half_edge( m, v0, v1 );
    if ( found == InvalidHalfEdgeID || !edge_info[found].is_spine( ))
    {
        return InvalidHalfEdgeID;
    }
    else
    {
        return found;
    }
}

    
}}
