//
//  geometric_properties.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "geometric_properties.h"
#include "polarize.h"
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <queue>
#include "Test.h"

using namespace CGLA;
using namespace std;
using namespace HMesh;
using namespace GLGraphics;
using namespace Procedural::Structure;

namespace Procedural{
    namespace Geometry{


Vec3d ring_barycenter ( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    vector< VertexID > vs;
    return ring_barycenter(m, h, vs);
}


Vec3d ring_barycenter ( Manifold& m, HalfEdgeID h, vector< VertexID > &vertices   )
{
    vector< HalfEdgeID > hes;
    return ring_barycenter(m , h, vertices, hes );
}


Vec3d ring_barycenter ( Manifold& m, HalfEdgeID h, vector< VertexID > &vertices,
                        vector< HalfEdgeID > &halfedges )
{
    assert(h != InvalidHalfEdgeID);
    
    Walker w = m.walker( h );
    Vec3d barycenter(0);
    
    while ( !w.full_circle( ))
    {
        vertices.push_back( w.vertex( ));
        halfedges.push_back(w.halfedge());
        barycenter += m.pos( w.vertex() );
        w = w.next().opp().next();
    }
    
    barycenter /= vertices.size();
    return barycenter;
}
        
void ring_vertices_and_halfedges ( Manifold& m, HalfEdgeID h, vector< HMesh::VertexID > &vertices,
                                    vector< HMesh::HalfEdgeID > &halfedges )
{
    ring_barycenter( m, h, vertices, halfedges );
}



int valence ( Manifold& m, VertexID v )
{
    Walker w = m.walker( v );
    for( ; !w.full_circle(); w = w.circulate_vertex_ccw());
    return w.no_steps();
}

bool is_singularity ( Manifold& m, VertexID v )
{
    return (! ( valence(m, v) == 4));
}
        

// this function is though to be used only for PAMs
Vec3d face_normal ( Manifold& m, FaceID f)
{
    assert( m.in_use( f ));
    
    Walker w = m.walker(f);
    vector< VertexID > verts;
    vector< Vec3d > pos;
    
    for ( ; !w.full_circle(); w = w.circulate_face_ccw( ))
    {
        verts.push_back( w.vertex( ));
        pos.push_back( m.pos( w.vertex( )));
    }
    
    // check if it is a quad or a triangle
    bool is_triangle = w.no_steps() == 3;
    if ( is_triangle )
    {
        auto n = cross( pos[1] - pos[0], pos[2] - pos[0] );
        //        n.normalize();
        return n;
    }
    else
    {
        auto n1 = cross( pos[1] - pos[0], pos[2] - pos[0] );
        auto n2 = cross( pos[2] - pos[0], pos[3] - pos[0] );
        auto n = ( n1 + n2 )/2.0;
        //        n.normalize();
        return n;
    }
}

Vec3d vertex_normal ( Manifold& m, VertexID v)
{
    Walker w = m.walker(v);
    Vec3d vn(0);
    for ( ; !w.full_circle(); w = w.circulate_vertex_ccw( ))
    {
        vn += face_normal(m, w.face());
    }
    
    vn / w.no_steps();
    //    vn.normalize();
    return vn;
}
        
bool is_2_neighbor_of_pole ( HMesh::Manifold& m, HMesh::VertexID v )
{
    bool itis = false;
    
    Walker w = m.walker( v );
    // check all of the vertices in its 1-ring
    for( ; !w.full_circle() && !itis; w = w.circulate_vertex_ccw( ))
    {
        itis = is_pole( m, w.vertex());
        if( itis ) continue;
        // circulate the 1-ring of the pointed vertex
        Walker wn = m.walker(w.vertex());
        for( ; !wn.full_circle() && !itis; wn = wn.circulate_vertex_ccw( ))
        {
            itis = is_pole( m, wn.vertex());
        }
    }
    return itis;
}
   
// should be in a separate class
double angle ( CGLA::Vec3d l, CGLA::Vec3d r)
{
    l.normalize();
    r.normalize();
    double cos = CGLA::dot( l, r );
    return acos(cos);
}
        
void vertex_distance_from_poles ( Manifold& m, VertexID v,
                                               HalfEdgeAttributeVector<EdgeInfo> edge_info,
                                               VertexAttributeVector<DistanceMetrics> &distances )
{
    // AS START, USE ONLY THE HOP METRIC
    int max_dist = numeric_limits<int>::min();
    std::queue< VertexID > poles;
    std::map< VertexID, DistanceMetrics > ds;
    // put all the poles into a pole queue
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        if( is_pole( m, *vit )) poles.push( *vit );
    }
    // while the queue is not empty
    while( !poles.empty( ))
    {
    //  set current_distance to 0
        int current_distance = 0;
    //  extract a pole from the queue.
        VertexID curr = poles.front();
        poles.pop();
        //  set 0 the distance
        ds[curr] = make_pair( current_distance, 0 );
        //  current_distance++
        current_distance++;
    //  take the pole 1 ring - and set the distance to current_distance.
        HalfEdgeID  he_helper       = m.walker( curr ).halfedge();
        VertexID    pole_neighbor   = m.walker( curr ).vertex();
        //  while the representative vertex of the next ring is not a pole
        int count = 0;
        while ( !is_pole( m, pole_neighbor ))
        {
//            cout << " looping for time : " << ( count++ ) << endl;
            HalfEdgeID ring_he = get_first_rib_edge( m, pole_neighbor, edge_info, true );
            vector< VertexID >      ring_vs;
            vector< HalfEdgeID >    ring_hes;
            ring_vertices_and_halfedges( m, ring_he, ring_vs, ring_hes );
            for( auto v : ring_vs )
            {
                //      if the distance for that ring of vertex is not set || is smaller than current_distance
                //          set the distance to current_distance
                if( ds.count(v) == 0 )
                {
                    ds[v] = make_pair( current_distance, 0 );
                }
                if( ds.count(v) >  0 )
                {
                    if(  ds[v].first > current_distance )
                    {
                        ds[v] = make_pair( current_distance, 0 );
                    }
                }
            }
            current_distance++;
    //      move to the next ring
            he_helper       = m.walker( he_helper ).next().opp().next().halfedge();
            pole_neighbor   = m.walker( he_helper ).vertex();
            
        }
    }
    cout << "finito";
    // debug
    for( auto fid : m.faces())
    {
        DebugRenderer::face_colors[fid] = Vec3f( 1.0, 1.0, 1.0 );
    }
    // vertices in some strange configurations are not reacheable with this procedure
    // fix that starting from those vertices, follow each of their edge chain
    // until the circle is clsoed or a pole is found.
    // a pole should be always found.
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        // take each vertex that does not have distance set
        if ( !(ds.count( *vit ) > 0 ))
        {
            Walker w = m.walker(*vit);
            for( ; !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                if( edge_info[w.halfedge()].is_rib( )) continue;
                Walker lw = m.walker(w.halfedge());
                int current_distance = 0;
                while ( !lw.full_circle( ) && !is_pole(m, lw.vertex()))
                {
                    lw = lw.next().opp().next();
                    current_distance++;
                }
                if ( is_pole(m, lw.vertex( )))
                {
                    if( !(ds.count( *vit ) > 0 ))
                       ds[*vit] = DistanceMetrics( current_distance, 0.0);
                    else
                    {
                        if( ds[*vit].first > current_distance )
                        {
                            ds[*vit].first = current_distance;
                        }
                    }
                }
            }
        }
        if( ds[*vit].first > max_dist ) max_dist = ds[*vit].first;
    }
    // debug
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        if ( !(ds.count( *vit ) > 0 ))
        {
            Walker w = m.walker(*vit);
            for(; !w.full_circle(); w = w.circulate_vertex_ccw())
                DebugRenderer::face_colors[w.face()] = Vec3f( 1.0, 1.0, 0.0);
        }
    }
    //
    
    
    
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        assert( m.in_use( *vit ));
        assert( ds.count( *vit ) > 0 );
        
        if( is_pole(m, *vit )) assert( ds[*vit].first == 0 );
        
        
        Vec3f color = color_ramp(ds[*vit].first, max_dist);
        DebugRenderer::vertex_colors[*vit] = color;
        Walker w = m.walker(*vit);
        for(; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            if( edge_info[w.halfedge()].is_rib() )
            {
                DebugRenderer::edge_colors[w.halfedge()]        = color;
                DebugRenderer::edge_colors[w.opp().halfedge()]  = color;
            }
        }

    }
}



}}