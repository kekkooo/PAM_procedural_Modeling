//
//  structural_helpers.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "structural_helpers.h"

#include "geometric_properties.h"

using namespace HMesh;
using namespace std;

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
        
HalfEdgeID get_first_rib_edge( Manifold& m, VertexID v0, HalfEdgeAttributeVector<EdgeInfo> edge_info, bool allow_junction )
{
    Walker      w   = m.walker( v0 );
    HalfEdgeID  he  = InvalidHalfEdgeID;
    for( ; !w.full_circle() && he == InvalidHalfEdgeID; w = w.circulate_vertex_ccw( ))
    {
        if( edge_info[w.halfedge()].is_rib() )
        {
            if( edge_info[w.halfedge()].is_junction() )
            {
                if( allow_junction ) { he = w.halfedge(); }
            }
            else
            {
                he = w.halfedge();
            }
        }
    }
    return he;
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
        

HalfEdgeID get_ring_vertices ( Manifold& m, VertexID v0, HalfEdgeAttributeVector<EdgeInfo> edge_info,
                         vector< HMesh::VertexID > &vertices, bool allow_junction)
{
    HalfEdgeID he = get_first_rib_edge( m, v0, edge_info, allow_junction );
    vertices.clear();
    Walker w = m.walker( he );
    for( ; !w.full_circle(); w = w.next().opp().next() )
    {
        vertices.push_back( w.vertex() );
    }
    return he;

}
        
void get_rings_from_pole ( HMesh::Manifold& m, HMesh::VertexID pole,
                           HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info,
                           std::vector<HMesh::HalfEdgeID> &rings,
                           int max_rings )
{
    if( !is_pole( m, pole )) return;
    rings.clear();
    int max_iter = max_rings < 0 ? numeric_limits<int>::max() : max_rings;
    // get the starters
    Walker w = m.walker( pole );
    while( !edge_info[w.next().halfedge()].is_junction()
        && !is_pole( m, w.vertex( )) && rings.size() < max_iter )
    {
        rings.push_back( w.next().halfedge() );
        w = w.next().opp().next();
    }

    rings.push_back( w.next().halfedge() );
    
}


        
        
RegionBoundaries FollowSpines ( Manifold& m, VertexID vid, HalfEdgeAttributeVector<EdgeInfo> edge_info )
{
    RegionBoundaries result = RegionBoundaries();
    
    if( is_pole( m, vid )) return result;

    HalfEdgeID  spine1 = InvalidHalfEdgeID,
                spine2 = InvalidHalfEdgeID;
    // find the two spine edge starters
    Walker w = m.walker( vid );
    for( ; !w.full_circle(); w = w.circulate_vertex_ccw())
    {
        if( edge_info[w.halfedge()].is_spine() )
        {
            if( spine1 == InvalidHalfEdgeID) { spine1 = w.halfedge(); }
            else spine2 = w.halfedge();
        }
    }
    // follow each of the edge paths until you reach a pole or a junction
    w = m.walker(spine1);
    while( !edge_info[w.next().halfedge()].is_junction() && !is_pole(m, w.vertex()) )
    {
        w = w.next().opp().next();
    }
    bool spine1_ends_to_junction = !edge_info[w.next().halfedge()].is_junction();
    // second starter
    w = m.walker(spine2);
    while( !edge_info[w.next().halfedge()].is_junction() && !is_pole(m, w.vertex()) )
    {
        w = w.next().opp().next();
    }
    bool spine2_ends_to_junction = !edge_info[w.next().halfedge()].is_junction();

    // return the case.
    result.FirstSpineEdge   = spine1;
    result.FirstEndType     = ( spine1_ends_to_junction ? EndType::Junction : EndType::Pole );
    result.FirstSpineEdge   = spine2;
    result.FirstEndType     = ( spine2_ends_to_junction ? EndType::Junction : EndType::Pole );

    return result;
}
        
// ALERT : IS NOT THAT EASY. YOU NEED TO BUILD THE GRAPH IN ORDER TO FIND LOOPS
HalfEdgeID BranchCutDirection( Manifold& m, VertexID vid, HalfEdgeAttributeVector<EdgeInfo> edge_info )
{
    HalfEdgeID result = InvalidHalfEdgeID;
    RegionBoundaries boundaries = FollowSpines( m, vid, edge_info );
    
    if( boundaries.FirstEndType == EndType::Pole ) result = boundaries.FirstSpineEdge;
    else if( boundaries.SecondEndType == EndType::Pole ) result = boundaries.SecondSpineEdge;
    
    return result;
}
        
void LabelJunctions ( Manifold& m, HalfEdgeAttributeVector<EdgeInfo> &edge_info )
{
    for(auto h: m.halfedges())
    {
        if(edge_info[h].is_rib())
        {
            number_rib_edges(m, edge_info, h);
            break;
        }
    }
}
        
/// gets the pole that is on the other end of the chain of opposite rib edges
VertexID get_other_end_pole  ( Manifold &m, VertexID pole, FaceID f)
{
    assert( pole != InvalidVertexID );
    assert( f    != InvalidFaceID   );
    assert( m.walker(f).vertex() == pole || m.walker(f).next().vertex() == pole ||
            m.walker(f).prev().vertex() == pole );
    assert( is_pole( m, pole ));
    
    Walker w = m.walker( f );
    for( ; w.vertex() != pole; w = w.next() );
    Walker right = w, left = w.next();
    while( !is_pole( m, left.vertex() ))
    {
        left    = left.next().opp().next();
        right   = right.prev().opp().prev();
    }
    assert( left.face() == right.face() );
    assert( left.next().halfedge() == right.halfedge());
    assert( is_pole(m, left.vertex( )));
    
    return left.vertex();
}
        
bool safe_to_add_branch( Manifold &m, VertexID v, size_t size, HalfEdgeAttributeVector<EdgeInfo> &edge_info ){
    assert( size > 0 );
    if( is_pole( m, v )) { return false; }
    if( Procedural::Geometry::is_neighbor_of_pole( m, v )) { return false; }
    if( valency( m, v ) != 4 ) { return false; }
    
    Walker w = m.walker( v );
    for( ; !w.full_circle() && !edge_info[w.halfedge()].is_rib(); w = w.circulate_vertex_ccw());
    assert(!w.full_circle()); // if true, probably edge_info is not initialized
    Walker left     = w;
    Walker right    = w.circulate_vertex_ccw().circulate_vertex_ccw();
    
    bool found_pole = false;
    for( size_t i = 0; i < size && !found_pole; ++i ){
        found_pole = found_pole || is_pole(m, left.vertex( )) || is_pole(m, right.vertex( ));
        left  = left.next().opp().next();
        right = right.next().opp().next();
    }
    return  ( !found_pole );
}

    
}}
