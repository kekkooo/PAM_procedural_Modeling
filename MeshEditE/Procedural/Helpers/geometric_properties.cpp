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
#include "math.h"

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
        
double dihedral_angle ( Vec3d a, Vec3d b, Vec3d c, Vec3d d )
{
    Vec3d u = b - a;     u.normalize();
    Vec3d v = c - a;     v.normalize();
    Vec3d w = d - a;     w.normalize();
    // classic formulation
//    Vec3d  uv = cross( u, v );
//    Vec3d  uw = cross( u, w );
//    double cs  = dot( uv, uw );
//    double ag  = acos( cs );
    
    // limit formulation
    double limit_cs = ( dot( u, u ) * dot( v, w )) - ( dot( u, w ) * dot( u, v ));
    double limit_ag = acos( limit_cs );
    
    // debug
    if( isnan( limit_cs ))
    {
        cout << a << endl;
        cout << b << endl;
        cout << c << endl;
        cout << d << endl;
    }
    
    assert( !isnan( limit_cs ));
    //

    
//    cout << " standard formulation " << cs << " # " << ag << endl
//         << " limit formulation "    << limit_cs << " # " << limit_ag << endl;
    
    return limit_cs;
}
        
void distance_from_poles ( Manifold& m, const HalfEdgeAttributeVector<EdgeInfo> edge_info,
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
        while ( !is_pole( m, pole_neighbor ))
        {
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
        distances[*vit] = ds[*vit];
    }
    
    // debug - useful
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
        //
    }
    cout << "finito";
}
        
        
void distance_from_junctions ( Manifold& m, const HalfEdgeAttributeVector<EdgeInfo> edge_info,
                                 VertexAttributeVector<DistanceMetrics> &distances )
{
    // AS START, USE ONLY THE HOP METRIC
    int max_dist = numeric_limits<int>::min();
    std::queue< VertexID > poles;
    std::map< VertexID, DistanceMetrics > ds;
    bool has_junctions = false;
    
    for( auto heid : m.halfedges())
    {
        int current_distance = 0;
        // skip the non junction edges
        if( !(edge_info[heid].is_junction( ))) continue;
        has_junctions = true;
        
        Walker w = m.walker( heid );
        if ( ds.count(w.vertex()) == 0 )
            ds[w.vertex()] = DistanceMetrics( current_distance, 0.0 );
        
        do{
            current_distance++;
            auto vid = w.next().vertex();
            if ( ds.count( vid ) == 0 )
                ds[ vid ] = DistanceMetrics( current_distance, 0.0 );
            else
            {
                if( ds[vid].first > current_distance )
                    ds[vid].first = current_distance;
            }
            w = w.next().next().opp();
        }while( !is_pole(m, w.next().vertex())
             && !edge_info[ w.next().next().halfedge()].is_junction() );
        
        if( is_pole(m, w.next().vertex( )))
            ds[ w.next().vertex( ) ] = DistanceMetrics( current_distance, 0.0 );
    }
    
    // all the following stuff is meaningful only if the mesh has junctions
    if(has_junctions)
    {
        // debug
        for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
        {
            if ( !(ds.count( *vit ) > 0 ))
            {
                Walker w = m.walker(*vit);
                for(; !w.full_circle(); w = w.circulate_vertex_ccw())
                    DebugRenderer::face_colors[w.face()] = Vec3f( 1.0, 1.0, 0.0);
            }
            else
            {
                max_dist = max_dist < ds[*vit].first ? ds[*vit].first : max_dist;
            }

        }
        //

        Vec3f blue( 0.0, 0.0, 1.0 );

        for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
        {
            if ( *vit == InvalidVertexID ) continue;
            assert( m.in_use( *vit ));
            assert( ds.count( *vit ) > 0 );
            
            Vec3f color = color_ramp(ds[*vit].first, max_dist);
            DebugRenderer::vertex_colors[*vit] = color;
            Walker w = m.walker(*vit);
            for(; !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                if( edge_info[w.halfedge()].is_rib() && !edge_info[w.halfedge()].is_junction() )
                {
                    DebugRenderer::edge_colors[w.halfedge()]        = color;
                    DebugRenderer::edge_colors[w.opp().halfedge()]  = color;
                }
                if( edge_info[w.halfedge()].is_junction() )
                {
                    DebugRenderer::edge_colors[w.halfedge()]        = blue;
                    DebugRenderer::edge_colors[w.opp().halfedge()]  = blue;
                }
            }
            
        }
    }
}
        
void distance_from_poles_and_junctions ( Manifold& m,
                                         const HalfEdgeAttributeVector<EdgeInfo> edge_info,
                                         VertexAttributeVector<DistanceMetrics> &distances )
{
    std::map< VertexID, DistanceMetrics > ds;
    
    // get distances from junctions
    //
    for( auto heid : m.halfedges())
    {
        int current_distance = 0;
        // skip the non junction edges
        if( !(edge_info[heid].is_junction( ))) continue;
        
        Walker w = m.walker( heid );
        if ( ds.count(w.vertex()) == 0 )
            ds[w.vertex()] = DistanceMetrics( current_distance, 0.0 );
        
        do{
            current_distance++;
            auto vid = w.next().vertex();
            if ( ds.count( vid ) == 0 )
                ds[ vid ] = DistanceMetrics( current_distance, 0.0 );
            else
            {
                if( ds[vid].first > current_distance )
                    ds[vid].first = current_distance;
            }
            w = w.next().next().opp();
        }while( !is_pole(m, w.next().vertex())
               && !edge_info[ w.next().next().halfedge()].is_junction() );
        
        if( is_pole(m, w.next().vertex( )))
            ds[ w.next().vertex( ) ] = DistanceMetrics( current_distance, 0.0 );
    }
    
    // add distances from poles
    //
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        if( !is_pole( m, *vit )) continue;
        VertexID curr = *vit;
        
        //  set current_distance to 0
        int current_distance = 0;
        
        //  set 0 the distance
        ds[curr] = make_pair( current_distance, 0 );
        //  current_distance++
        current_distance++;
        //  take the pole 1 ring - and set the distance to current_distance.
        HalfEdgeID  he_helper       = m.walker( curr ).halfedge();
        VertexID    pole_neighbor   = m.walker( curr ).vertex();
        //  while the representative vertex of the next ring is not a pole
        while ( !is_pole( m, pole_neighbor ))
        {
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
    int max_dist = numeric_limits<int>::min();
    
    // debug
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        if ( !(ds.count( *vit ) > 0 ))
        {
            Walker w = m.walker(*vit);
            for(; !w.full_circle(); w = w.circulate_vertex_ccw())
                DebugRenderer::face_colors[w.face()] = Vec3f( 1.0, 1.0, 0.0);
        }
        else
        {
            max_dist = max_dist < ds[*vit].first ? ds[*vit].first : max_dist;
        }
        distances[*vit] = ds[*vit];
        
    }
    
    Vec3f blue( 0.0, 0.0, 1.0 );
    
    for( auto vit = m.vertices().begin(); vit != m.vertices().end(); ++vit )
    {
        assert( m.in_use( *vit ));
        assert( ds.count( *vit ) > 0 );
        
        Vec3f color = color_ramp(ds[*vit].first, max_dist);
        DebugRenderer::vertex_colors[*vit] = color;
        Walker w = m.walker(*vit);
        for(; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            if( edge_info[w.halfedge()].is_rib() && !edge_info[w.halfedge()].is_junction() )
            {
                DebugRenderer::edge_colors[w.halfedge()]        = color;
                DebugRenderer::edge_colors[w.opp().halfedge()]  = color;
            }
            if( edge_info[w.halfedge()].is_junction() )
            {
                DebugRenderer::edge_colors[w.halfedge()]        = blue;
                DebugRenderer::edge_colors[w.opp().halfedge()]  = blue;
            }
        }
    }
}
        
Vec3f get_angle_color( double coseno )
{
    return color_ramp2( (int)( ( 1.0 + coseno ) * 20.0 ) , 20 );
}
        
void dihedral_angles ( Manifold& m, HalfEdgeAttributeVector<EdgeInfo> edge_info,
                       VertexAttributeVector<double> &angles, EdgeType edge_type )
{
    for( auto a : m.vertices() )
    {
        // skip poles.
        if( is_pole( m, a )) continue;
        
        // select the correct set of vertices depending on the edge type
        VertexID    b1,  b2, c, d;
        HalfEdgeID  he1, he2;
        Walker      w           = m.walker(a);
        bool        rib_start   = edge_info[w.halfedge()].is_rib();

        VertexID north = rib_start ? w.circulate_vertex_ccw().vertex() : w.vertex(),
                 east  = rib_start ? w.vertex() : w.circulate_vertex_cw().vertex(),
                 south = rib_start ? w.circulate_vertex_cw().vertex()
                                   : w.circulate_vertex_ccw().circulate_vertex_ccw().vertex(),
                 west  = rib_start ? w.circulate_vertex_ccw().circulate_vertex_ccw().vertex()
                                   : w.circulate_vertex_ccw().vertex();
        
        if( edge_type == RIB || edge_type == RIB_JUNCTION )
        {
            b1  = east;  b2 = west;
            c   = north; d  = south;
        }
        else
        {
            b1 = north;  b2 = south;
            c  = west;   d  = east;
        }
        he1 = find_half_edge( m, a, b1 );
        he2 = find_half_edge( m, a, b2 );
        
        assert( m.in_use( a  ));
        assert( m.in_use( b1 ));
        assert( m.in_use( b2 ));
        assert( m.in_use( c  ));
        assert( m.in_use( d  ));
        
        double cos1  = dihedral_angle( m.pos( a ), m.pos( b1 ), m.pos( c ), m.pos( d ));
        double cos2  = dihedral_angle( m.pos( a ), m.pos( b2 ), m.pos( c ), m.pos( d ));
        double cosv  = ( cos1 + cos2 ) / 2.0;
        
        angles[a] = cosv;
        
        DebugRenderer::edge_colors[he1]                             = get_angle_color( cos1 );
        DebugRenderer::edge_colors[m.walker(he1).opp().halfedge()]  = get_angle_color( cos1 );
        DebugRenderer::edge_colors[he2]                             = get_angle_color( cos2 );
        DebugRenderer::edge_colors[m.walker(he2).opp().halfedge()]  = get_angle_color( cos2 );
        DebugRenderer::vertex_colors[a]                             = get_angle_color( cosv );
    }
}
        
void dihedral_angles ( Manifold& m, const HalfEdgeAttributeVector<EdgeInfo> edge_info,
                       VertexAttributeVector<double> &angles )
{
    VertexAttributeVector<double> rib_angles, spine_angles;
    dihedral_angles( m, edge_info, rib_angles, RIB );
    dihedral_angles( m, edge_info, spine_angles, SPINE );
    for( auto vid : m.vertices( ))
    {
        angles[vid] = ( rib_angles[vid] + spine_angles[vid] ) / 2.0;
        DebugRenderer::vertex_colors[vid] = get_angle_color( angles[vid] );
    }
    
}
        
double edge_length ( Manifold& m,
                     HalfEdgeAttributeVector<double> lengths )
{
    std::map< HalfEdgeID, bool > visited;
    double                       length_sum = 0.0;
    unsigned int                 count      = 0;
    for( HalfEdgeID he : m.halfedges() )
    {
        if( !visited[he] )
        {
            Walker w = m.walker(he);
            visited[he]                 = true;
            visited[w.opp().halfedge()] = true;
            double l = ( m.pos( w.vertex( )) - m.pos( w.prev().vertex( ))).length();
            lengths[he]                 = l;
            lengths[w.opp().halfedge()] = l;
            length_sum                 += l;
            ++count;
        }
    }
    return length_sum / count;    
}
        
        
double mean_length_of_outoing_he ( Manifold& m, VertexID vertex )
{
    double length   = 0.0;
    Vec3d  curr     = m.pos( vertex );
    Walker w        = m.walker( vertex );
    for ( ; !w.full_circle(); w = w.circulate_vertex_ccw( ))
    {
        length += ( curr - m.pos( w.vertex( ))).length();
    }
    length /= w.no_steps();
    return  length;
}




}}