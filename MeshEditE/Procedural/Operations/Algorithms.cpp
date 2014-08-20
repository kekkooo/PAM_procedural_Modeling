//
//  Algorithms.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "Algorithms.h"

using namespace std;
using namespace CGLA;
using namespace HMesh;

namespace Procedural{
    namespace Operations{
        namespace Algorithms    {
    
void selected_vertices_inverse_distance_laplacian( HMesh::Manifold&m, vector< VertexID > vs )
{
    map< VertexID, Vec3d > new_pos;
    
    for( VertexID vid : vs )
    {
        double weight_sum   = 0.0;
        Walker w            = m.walker( vid );
        Vec3d vid_pos       = m.pos( vid );
        new_pos[vid]        = Vec3d( 0 );
        
        // iterate the one-ring of the vertex vid
        for( ; !w.full_circle( ); w = w.circulate_vertex_ccw( ) )
        {
            Vec3d   neighbor_pos    = m.pos( w.vertex() );
            double  length          = ( vid_pos - neighbor_pos ).length();
            
            assert( vid_pos != neighbor_pos );
            assert( !isnan( length ));
            
            double  weight   = 1.0 / length;
            new_pos[vid]    += ( neighbor_pos * weight );
            weight_sum      += weight;
        }
        new_pos[vid] /= weight_sum;
    }
    // update coordinates
    for( VertexID vid : vs )  { m.pos( vid ) = new_pos[vid]; }
}
            
void selected_vertices_cotangent_weights_laplacian( HMesh::Manifold&m, vector< VertexID > vs )
{
    map< VertexID, Vec3d > new_pos;
    
    for( VertexID vid : m.vertices( ) )
    {
        double weight_sum   = 0.0;
        Walker w            = m.walker( vid );
        Vec3d vid_pos       = m.pos(vid);
        new_pos[vid]        = Vec3d( 0 );
        
        for ( ; !w.full_circle(); w = w.circulate_vertex_ccw( ))
        {
            Vec3d   neighbor_pos    = m.pos(w.vertex());
            double  length          = ( vid_pos - neighbor_pos ).length();
            assert( !isnan( length ) && length > 0 );
            HalfEdgeID  curr_he     = w.halfedge();
            HalfEdgeID  twin_he     = w.opp().halfedge();
            
            Vec3d   face_1_centroid(0);
            Vec3d   face_2_centroid(0);
            
            Vec3d   edge_midpoint( ( vid_pos[0] + neighbor_pos[0] )/2.0,
                                   ( vid_pos[1] + neighbor_pos[1] )/2.0,
                                   ( vid_pos[2] + neighbor_pos[2] )/2.0);
            Walker wf1 = m.walker( curr_he );
            Walker wf2 = m.walker( twin_he );
            // calculate centroid of the first face
            for ( ; !wf1.full_circle(); wf1 = wf1.next( )) {
                face_1_centroid += m.pos( wf1.vertex( ));
            }
            face_1_centroid /= (double)wf1.no_steps();
            // calculate centroid of the second face
            for ( ; !wf2.full_circle(); wf2 = wf2.next( )) {
                face_2_centroid += m.pos( wf2.vertex( ));
            }
            face_2_centroid /= (double)wf2.no_steps();
            // calculate heights as distance between the edge midpoint and the face's centroid
            double h1 = ( edge_midpoint - face_1_centroid ).length();
            double h2 = ( edge_midpoint - face_2_centroid ).length();
            assert( !isnan( h1 ) && h1 > 0);
            assert( !isnan( h2 ) && h2 > 0);
            // calculate weight and new position of the vertex
            double weight    = ( h1 + h2 ) / length;
            weight_sum      += weight;
            new_pos[vid]    += ( neighbor_pos * weight );
        }
        
        new_pos[vid] /= weight_sum;
    }
    // update coordinates
    for( VertexID vid : vs )  { m.pos( vid ) = new_pos[vid]; }
}
            
            
// weighted laplacian smoothing with inverse distance as weight
void inverse_distance_laplacian_smoothing( HMesh::Manifold& am )
{
    map< VertexID, Vec3d > new_pos;
    
    for( VertexID vid : am.vertices( ) )
    {
        double weight_sum   = 0.0;
        Walker w            = am.walker( vid );
        Vec3d vid_pos       = am.pos( vid );
        new_pos[vid]        = Vec3d( 0 );
        
        // iterate the one-ring of the vertex vid
        for( ; !w.full_circle( ); w = w.circulate_vertex_ccw( ) )
        {
            Vec3d   neighbor_pos    = am.pos( w.vertex() );
            double  length          = ( vid_pos - neighbor_pos ).length();
            
            assert( vid_pos != neighbor_pos );
            assert( !isnan( length ));
            
            double  weight          = 1.0 / length;
            new_pos[vid]   += ( neighbor_pos * weight );
            weight_sum += weight;
        }
        
        new_pos[vid] /= weight_sum;
    }
    // update coordinates
    for( VertexID vid : am.vertices( ) )  { am.pos( vid ) = new_pos[vid]; }
}

// weighted laplacian smoothing with cotangets as weight
void cotangent_weights_laplacian_smoothing( HMesh::Manifold& am)
{
    map< VertexID, Vec3d > new_pos;
    
    for( VertexID vid : am.vertices( ) )
    {
        double weight_sum   = 0.0;
        Walker w            = am.walker( vid );
        Vec3d vid_pos       = am.pos(vid);
        new_pos[vid]        = Vec3d( 0 );
        
        for ( ; !w.full_circle(); w = w.circulate_vertex_ccw( ))
        {
            Vec3d   neighbor_pos    = am.pos(w.vertex());
            double  length          = ( vid_pos - neighbor_pos ).length();
            assert( !isnan( length ) && length > 0 );
            HalfEdgeID  curr_he     = w.halfedge();
            HalfEdgeID  twin_he     = w.opp().halfedge();
            
            Vec3d   face_1_centroid(0);
            Vec3d   face_2_centroid(0);
            
            Vec3d   edge_midpoint( ( vid_pos[0] + neighbor_pos[0] )/2.0,
                                  ( vid_pos[1] + neighbor_pos[1] )/2.0,
                                  ( vid_pos[2] + neighbor_pos[2] )/2.0);
            Walker wf1 = am.walker(curr_he);
            Walker wf2 = am.walker(twin_he);
            // calculate centroid of the first face
            for ( ; !wf1.full_circle(); wf1 = wf1.next( )) {
                face_1_centroid += am.pos( wf1.vertex( ));
            }
            face_1_centroid /= (double)wf1.no_steps();
            // calculate centroid of the second face
            for ( ; !wf2.full_circle(); wf2 = wf2.next( )) {
                face_2_centroid += am.pos( wf2.vertex( ));
            }
            face_2_centroid /= (double)wf2.no_steps();
            // calculate heights as distance between the edge midpoint and the face's centroid
            double h1 = ( edge_midpoint - face_1_centroid ).length();
            double h2 = ( edge_midpoint - face_2_centroid ).length();
            assert( !isnan( h1 ) && h1 > 0);
            assert( !isnan( h2 ) && h2 > 0);
            // calculate weight and new position of the vertex
            double weight    = ( h1 + h2 ) / length;
            weight_sum      += weight;
            new_pos[vid]    += ( neighbor_pos * weight );
        }
        
        new_pos[vid] /= weight_sum;
    }
    // update coordinates
    for( VertexID vid : am.vertices( ) )  { am.pos( vid ) = new_pos[vid]; }
}
            
            
void along_spines ( Manifold& m, HalfEdgeAttributeVector<EdgeInfo> &edge_info )
{
    map< VertexID, Vec3d > new_pos;
    
    for( VertexID vid : m.vertices( ) )
    {
        double weight_sum   = 0.0;
        Walker w            = m.walker( vid );
        Vec3d vid_pos       = m.pos( vid );
        new_pos[vid]        = Vec3d( 0 );
        
        // iterate the one-ring of the vertex vid
        for( ; !w.full_circle( ); w = w.circulate_vertex_ccw( ) )
        {
            // skip all ribs
            if( edge_info[ w.halfedge() ].is_rib()) continue;
            Vec3d   neighbor_pos    = m.pos( w.vertex() );
            double  length          = ( vid_pos - neighbor_pos ).length();
            
            assert( vid_pos != neighbor_pos );
            assert( !isnan( length ));
            
            double  weight          = 1.0 / length;
            new_pos[vid]   += ( neighbor_pos * weight );
            weight_sum += weight;
        }
        
        new_pos[vid] /= weight_sum;
    }
    // update coordinates
    for( VertexID vid : m.vertices( ) )  { m.pos( vid ) = new_pos[vid]; }
}

}}} // namespaces closing