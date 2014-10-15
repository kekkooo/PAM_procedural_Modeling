//
//  patch_mapping.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 06/10/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "patch_mapping.h"
#include "polarize.h"
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include <queue>
#include <unordered_set>

using namespace HMesh;
using namespace std;
using namespace Procedural::Geometry;
using namespace GLGraphics;
using namespace CGLA;


bool is_corner ( Manifold& m, VertexID v, HalfEdgeAttributeVector<int> &layout )
{
    Walker w = m.walker( v );
    vector<HalfEdgeID> out;
    for (; !w.full_circle(); w = w.circulate_vertex_ccw()) { out.push_back( w.halfedge()); }
    if( out.size() != 4 || is_pole( m, v ) )
    {
        return true;
    }
    else
    {
        // if all outgoing edges are part of the layout it is a corner
        bool is_corner = true;
        for( HalfEdgeID heid : out ) { is_corner = is_corner && layout[heid]; }
        return is_corner;
    }
}

void build_patches( Manifold& m, FaceAttributeVector<int> &face_to_patch )
{
    HalfEdgeAttributeVector<int> is_layout, dont_start_from_here;
    map< HalfEdgeID, int > patch_boundaries;
    for( FaceID fid      : m.faces( ))      { face_to_patch[fid]    = -1; }
    for( HalfEdgeID heid : m.halfedges( ))  { is_layout[heid]       = 0;
                                              dont_start_from_here[heid] = 0; }
    
    std::vector<VertexID> singularities;
    std::vector<VertexID>    corners;
    
    // 1) collect all singularities that are not poles
    for( VertexID vid : m.vertices( ))
    {
        if( is_pole( m, vid )) continue;
        if( is_singularity( m, vid )) { singularities.push_back( vid ); }
    }
    // 2) for all singularities
    for( VertexID s : singularities )
    {
    // for each outgoing edge
        Walker w = m.walker( s );
        for (; !w.full_circle(); w = w.circulate_vertex_ccw( ))
        {
            assert( is_layout[w.halfedge()] == is_layout[w.opp().halfedge()] );
            // skip if the outgoing edge is marked as layout ( means it was visited )
            if( is_layout[w.halfedge()] ) continue;
            Walker w_chain = m.walker( w.halfedge( ));
            // follow the edge until it reaches a singularity
            while( !is_singularity( m, w_chain.vertex()) )
            {
                //  at each step set the current halfedge and its twin as part of layout
                is_layout[ w_chain.halfedge() ]       = 1;
                is_layout[ w_chain.opp().halfedge() ] = 1;
                w_chain = w_chain.next().opp().next();
            }
            is_layout[ w_chain.halfedge() ]       = 1;
            is_layout[ w_chain.opp().halfedge() ] = 1;
        }
    }
    
    for( HalfEdgeID he : m.halfedges() )
    {
        assert( is_layout[ he ] == is_layout[m.walker(he).opp().halfedge()] );
        
    }
    
    
    // - FIND ALL THE CORNERS
    // put the poles as first in order to start from them
    // since turning at the poles is tricky
    for( VertexID vid : m.vertices())
    {
        if( is_pole( m, vid )) { corners.push_back( vid ); }
    }
    for( VertexID vid : m.vertices())
    {
        if( !is_pole(m, vid))
        {
            if( is_corner( m, vid, is_layout )) { corners.push_back( vid ); }
        }
    }
    
    
    // BEFORE GOING FURTHER CHECK USING DEBUG RENDERER
    for( auto heid : m.halfedges() )
    {
        if( is_layout[heid])
            DebugRenderer::edge_colors[heid] = Vec3f( 0.0, 0.0, 1.0);
        else
            DebugRenderer::edge_colors[heid] = Vec3f( 0.0, 0.0, 0.0);
    }
//    for( FaceID fid : m.faces() )
//    {
//        DebugRenderer::face_colors[fid] = Vec3f( 0.0, 0.0, 0.0);
//    }
    for( VertexID c : corners )
    {
        DebugRenderer::vertex_colors[c] = Vec3f( 1.0, 1.0, 1.0 );
    }

    int current_patch = 1;
    // - ASSIGN PATCH ID TO EACH BORDER FACE
    for( VertexID C : corners )
    {
        // skip the poles
//        if( is_pole( m, C )) continue;
        Walker w = m.walker( C );
        vector<HalfEdgeID> out;
        for (; !w.full_circle(); w = w.circulate_vertex_ccw()) { out.push_back( w.halfedge()); }
        // for each outgoing edge
        for( HalfEdgeID he : out )
        {
            bool last_vertex_was_corner = false;
            // skip if visited
            if( dont_start_from_here[he] == 1 ) continue;
            if( !is_layout[he])                 continue;
            // mark as visited
            dont_start_from_here[he] = 1;
            // set visited->face's patch as current_patch
            Walker hew = m.walker( he );
            assert( hew.halfedge() == he );
            bool done = false;
            //      iterate until C is reached
            while ( !done )
            {
                cout << "last vertex was corner?" << last_vertex_was_corner << endl;
                assert(is_layout[hew.halfedge()] == 1);
                DebugRenderer::edge_colors[hew.halfedge()] = Vec3f( 1.0, 0.0, 0.0);
                face_to_patch[ hew.face() ] = current_patch;
                done = hew.vertex() == C;
                // follow until a corner is reached
                if( is_corner( hew.vertex(), corners ))
//                if( is_corner( m, hew.vertex(), is_layout ))
                {
                    DebugRenderer::vertex_colors[hew.vertex()] = Vec3f( 0.0, 1.0, 0.0);
                    // turn left
                    hew = hew.next();
                    // mark outgoing edge as visited
                    dont_start_from_here[hew.halfedge()] = 1;
                    last_vertex_was_corner = true;
                }
                else
                {
                    hew = hew.next().opp().next();
                    last_vertex_was_corner = false;
                }
            }
            face_to_patch[ hew.face() ] = current_patch;
            ++current_patch;
        }
        
    }
    
    // SANITY CHECK - check if the two sides of a layout edge have different patch
    for( HalfEdgeID he : m.halfedges() )
    {
        if( is_layout[he] == 1 )
        {
            Walker w = m.walker( he );
            assert( face_to_patch[w.face()] != face_to_patch[w.opp().face()] );
        }
    }

    
    // - SPREAD THE IDs
    // collect all the faces with  patch id
    std::queue<FaceID> to_visit;
    for( FaceID f : m.faces())
    {
        if (face_to_patch[f] != -1 ) { to_visit.push(f); }
    }
    
    while( !to_visit.empty() )
    {
        FaceID f = to_visit.front();
        to_visit.pop();
        Walker w = m.walker(f);
        for (; !w.full_circle(); w = w.next())
        {
            assert( w.face() == f);
            if( face_to_patch[w.opp().face()] == -1 )
            {
                face_to_patch[w.opp().face()] = face_to_patch[ f ];
                to_visit.push( w.opp().face() );
            }
        }
    }
    
    // set all the triangles to the same patch
    for( VertexID corner : corners )
    {
        bool pole = is_pole(m, corner);
        if( pole )
        {
            Walker w = m.walker(corner);
            for( ; !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                face_to_patch[w.face()] = 0;
            }
        }
        
    }
    
    
    int no_poles = 0;
    int no_val3  = 0;
    int no_val5  = 0;
    int no_val6p = 0;
    int no_val6p_no_pole = 0;

    
    for( VertexID v : m.vertices() )
    {
        bool pole = is_pole(m, v);
        int val = valency(m, v);
        if( val == 3) ++no_val3;
        if( val == 5) ++no_val5;
        if( val > 5)  ++no_val6p;
        if( !pole )
        {
            if( val > 5 )
                ++no_val6p_no_pole;
        }
        if( pole ) { ++no_poles; }
    }
    
    
    std::unordered_set< int > patches;
    
    // draw with colors and save the number of patches
    for( FaceID fid : m.faces() )
    {
//        if( face_to_patch[fid] != -1 )
        DebugRenderer::face_colors[fid] = get_color(face_to_patch[fid]);
        patches.insert( face_to_patch[fid] );
    }

    cout << "#poles : " << no_poles << " #patches : " << patches.size() << endl ;
    cout << "#val 3  :          " << no_val3 << endl;
    cout << "#val 5  :          " << no_val5 << endl;
    cout << "#val 6+ :          " << no_val6p << endl;
    cout << "#val 6+ no pole :  " << no_val6p_no_pole << endl;

    
}