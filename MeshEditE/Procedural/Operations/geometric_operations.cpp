//
//  geometric_operations.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "geometric_operations.h"
#include "polarize.h"
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include "test.h"

using namespace HMesh;
using namespace CGLA;
using namespace std;
using namespace Procedural::Geometry;

namespace Procedural{
    namespace Operations{
        namespace Geometric{

            
void add_noise ( HMesh::Manifold& m, int vertex_type_flags, double ratio, double cutoff )
{
    for ( VertexID vid : m.vertices( ))
    {
        assert( m.in_use(vid));
        if( ( vertex_type_flags & VertexType::POLE ) && is_pole( m, vid ))
        {
            //            Vec3d dir           = simple_random_direction(m, vid) * ratio * 3.0;
            Vec3d dir           = alt_simple_random_direction( m, vid );
            Vec3d polar_normal  = vertex_normal( m, vid );
            
            if( dot( dir, polar_normal ) < 0 ) dir = -dir;
            
            if( cutoff > 0.0 )
            {
                if ( dir.length() > cutoff )
                {
                    dir.normalize();
                    dir *= cutoff;
                }
            }
            
            double scaling_ratio = (( rand() % 100 ) - 50.0 )/ 200.0 + 1.0;
            
            cout << "scaling " << scaling_ratio << endl;
            extrude_pole( m, vid, dir, true, scaling_ratio );
            
            //            move_vertex(m, vid, dir);
        }
        else if( vertex_type_flags & VertexType::REGULAR )
        {
            Vec3d dir = simple_random_direction( m, vid, 3.0 ) * ratio;
            if( cutoff > 0.0 )
            {
                if ( dir.length() > cutoff )
                {
                    dir.normalize();
                    dir *= cutoff;
                }
            }
            
            move_vertex( m, vid, dir );
        }
        //        m.cleanup();
    }
    
}


// h must be a spine edge. Splits the specified half edge and all of it "parallel" edges,
// that are all the edges reachable as current_he.next().next().opp()
vector< VertexID > split_ring_of_quads( HMesh::Manifold& m, HMesh::HalfEdgeID h)
{
    vector< VertexID > inserted_vertices;
    // cannot split invalid half edge
    if( h == InvalidHalfEdgeID ) return inserted_vertices;
    
    vector< FaceID > faces_to_split;
    Walker w = m.walker(h);
    // split all the edges and keep track of all the faces that must be split
    do
    {
        inserted_vertices.push_back( m.split_edge(w.halfedge()) );
        faces_to_split.push_back(w.face());
        w = w.next().next().opp();
    }while (w.halfedge() != h);
    
    // split the faces
    for (int i = 0; i < inserted_vertices.size() - 1; i++)
    {
        m.split_face_by_edge(faces_to_split[i], inserted_vertices[i], inserted_vertices[i+1]);
    }
    // split the last face ( the one that closes the loop
    m.split_face_by_edge(faces_to_split.back(), inserted_vertices[inserted_vertices.size()-1], inserted_vertices.front());
    
    return inserted_vertices;
}

// h must be a rib edge. splits each rib edge that lay in the same path "pole-to-pole" of h
vector< VertexID > split_from_pole_to_pole( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    vector< VertexID > vertices;
    // cannot split invalid half edge
    if( h == InvalidHalfEdgeID ) return vertices;
    
    vector< FaceID >        faces_to_split;
    vector< HalfEdgeID >    edges_to_split;
    
    Walker w = m.walker(h);
    
    // reach one of the poles
    while( !is_pole(m, w.next().vertex()) )
    {
        w = w.next().next().opp();
    }
    assert( is_pole(m, w.next().vertex( )));
    
    vertices.push_back( w.next().vertex( ));
    faces_to_split.push_back( w.face( ));
    
    // change direction in order to move far from the pole the Walker has reached
    w = w.opp();
    
    // loop until the walker reaches the opposite pole. Saving edges and faces to split
    while( !is_pole(m, w.next().vertex()) )
    {
        edges_to_split.push_back(w.halfedge());
        faces_to_split.push_back( w.face());
        w = w.next().next().opp();
    }
    // add the last edge and face
    edges_to_split.push_back(w.halfedge());
    faces_to_split.push_back(w.face());
    
    // split all the edges keeping track of the vertices added to the mesh
    for( auto he_id : edges_to_split)
    {
        vertices.push_back( m.split_edge( he_id ));
    }
    // add the last vertex ( the pole )
    vertices.push_back(w.next().vertex());
    
    // create the new faces accordin to the new vertices
    for (int i = 0; i < vertices.size() - 1; i++)
    {
        m.split_face_by_edge(faces_to_split[i], vertices[i], vertices[i+1]);
    }
    
    return vertices;
}

void move_vertex ( HMesh::Manifold& m, HMesh::VertexID v, Vec3d vector_to_apply )
{
    m.pos( v ) += vector_to_apply;
}


// extrudes one vertex. It works differently if the vertex is a pole or not
// if length is NaN then the length is calculate as the average outgoing edge's length
void extrude_pole ( HMesh::Manifold& m, HMesh::VertexID v, double length,
                   bool pole_add_quads, double radius_preserving_ratio )
{
    //  early ternination
    if( !is_pole( m, v )) return;
    
    // if length is NaN calculate the average outgoing edge's length
    if( isnan(length) )
    {
        length = 0.0;
        Vec3d curr  = m.pos(v);
        Walker w    = m.walker(v);
        for (; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            length += ( curr - m.pos(w.vertex())).length();
        }
        length /= w.no_steps();
    }
    // move the vertex and store the list of the vertices that are added during the split
    // of the triangles fan
    Vec3d vn    = vertex_normal(m, v);
    vn.normalize();
    move_vertex( m, v, length * vn );
    
    // if the vertex is a pole and the "add_quads" option is selected it extrudes adding quads
    if( pole_add_quads )
    {
        add_ring_around_pole( m, v, radius_preserving_ratio );
    }
}

void extrude_pole( HMesh::Manifold&m, HMesh::VertexID v, Vec3d direction, bool add_quads, double scaling )
{
    assert( is_pole( m, v ));
    move_vertex( m, v, direction );

    if( add_quads )
    {
        add_ring_around_pole( m, v, scaling );
    }
}


void scale_ring_radius ( HMesh::Manifold& m, HMesh::HalfEdgeID h, double ratio )
{
    vector< VertexID > vertices;
    Vec3d barycenter = ring_barycenter( m, h, vertices );
    for( auto v : vertices)
    {
        m.pos(v) = barycenter + ( m.pos( v ) - barycenter ) * ratio;
    }
}


void set_ring_radius ( Manifold& m, HalfEdgeID h, double radius )
{
    vector< VertexID > vertices(0);
    Vec3d barycenter = ring_barycenter( m, h, vertices );
    for( auto v : vertices)
    {
        Vec3d dir = ( m.pos( v ) - barycenter );
        dir.normalize();
        m.pos( v ) = barycenter + ( dir * radius );
    }
}

// the radii vector must contain at radii[0] the radius that has to be assigned to the ring at the vertex pointed to by h
void set_ring_radius ( Manifold& m, HalfEdgeID h, vector<double> radii )
{
    vector< VertexID > vertices;
    Vec3d barycenter = ring_barycenter( m, h, vertices );
    DistanceDirectionVector ddv;
    for( auto v : vertices)
    {
        Vec3d dir = ( m.pos( v ) - barycenter );
        double length = dir.length();
        dir.normalize();
        
        ddv.push_back( make_tuple( v, length, dir ));
    }
    
    set_ring_radius( m, barycenter, ddv );
}


void set_ring_radius ( HMesh::Manifold& m, Vec3d barycenter,  DistanceDirectionVector ddv )
{
    for( auto triple : ddv )
    {
        cout << " dist " << get<1>( triple ) << endl;
        m.pos( get<0>( triple )) = barycenter + ( get<2>( triple ) * get<1>( triple ));
    }
}


void flatten_pole ( Manifold& m, VertexID pole )
{
    if( !is_pole(m, pole)) return;
    
    // take all the neighbors,
    // find the combination of three that creates the biggest triangle
    // take those 3 points as plane
    // move the pole in its projection onto the plane.
    assert(false);
}
            
void add_ring_around_pole ( Manifold& m, VertexID pole, double scaling )
{
    // Add the other ring of rib edge loops
    HMesh::VertexAttributeVector<int> sel;
    for( VertexID vid : m.vertices() ) { sel[vid] = 0; }
    sel[pole] = 1;
    refine_poles( m, sel );
    
    // find the starting half_edges, the rings, their barycenters, and the associated radii
    vector< VertexID >  sourceRing,     destRing;
    vector<double>      sourceRadii,    destRadii;
    Walker w = m.walker(pole);
    // starters
    HalfEdgeID  destStart       = w.next().halfedge();
    HalfEdgeID  sourceStart     = w.next().opp().next().next().halfedge();
    // barycenters
    Vec3d       destBarycenter  = ring_barycenter( m, destStart,    destRing );
    Vec3d       sourceBarycenter= ring_barycenter( m, sourceStart,  sourceRing );
    
    // TEST TEST TEST
    assert( destRing.size() == sourceRing.size( ));
    for( int j = 0; j < destRing.size(); ++j )
    {
        bool found = false;
        Walker www = m.walker( destRing[j] );
        for ( ; !www.full_circle() && ! found; www = www.circulate_vertex_ccw( ))
        {
            found = www.vertex() == sourceRing[ j ];
        }
        assert(found);
    }
    // #############
    
    for( int j = 0; j < destRing.size(); ++j )
    {
        Vec3d   dir     = m.pos( destRing[j] )   - destBarycenter;
        dir.normalize();
        double  length  = ( m.pos( sourceRing[j] ) - sourceBarycenter ).length();
        m.pos( destRing[j] ) = destBarycenter + dir * length * scaling;
    }
}
    

}}}
