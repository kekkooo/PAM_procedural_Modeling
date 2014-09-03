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
#include <MeshEditE/Procedural/Helpers/plane.h>
#include <MeshEditE/Procedural/Operations/Algorithms.h>
#include <MeshEditE/Procedural/Operations/Algorithms.h>
#include "test.h"

using namespace HMesh;
using namespace CGLA;
using namespace std;
using namespace Procedural::Geometry;
using namespace Procedural::Operations::Algorithms;

namespace Procedural{
    namespace Operations{
        namespace Geometric{

            
void add_noise_to_vertex( Manifold& m, VertexID vid, int vertex_type_flags, double ratio, double cutoff )
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
        
//        cout << "scaling " << scaling_ratio << endl;
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
}


void add_noise ( HMesh::Manifold& m, int vertex_type_flags, double ratio, double cutoff,
                 vector< VertexID > vertices )
{
    for( auto vid : vertices )
    {
        add_noise_to_vertex(m, vid, vertex_type_flags, ratio, cutoff );
    }
}
    
            
void add_noise ( HMesh::Manifold& m, int vertex_type_flags, double ratio, double cutoff )
{
    for( auto vid : m.vertices( ))
    {
        add_noise_to_vertex(m, vid, vertex_type_flags, ratio, cutoff );
    }

}


// h must be a spine edge. Splits the specified half edge and all of it "parallel" edges,
// that are all the edges reachable as current_he.next().next().opp()
vector< VertexID > split_ring_of_quads( HMesh::Manifold& m, HMesh::HalfEdgeID h )
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
            
void split_ring_of_quads ( HMesh::Manifold& m, HMesh::HalfEdgeID h, int slices )
{
    typedef pair<Vec3d, Vec3d> point_dir;
    if( h == InvalidHalfEdgeID || slices <= 1 ) return;
    
    vector< VertexID >          inserted_vertices;
    vector< vector< Vec3d > >   inserted_vertices_position;
    vector< FaceID >            faces_to_split;
    vector< point_dir >         points_dirs;

    Walker w = m.walker(h);
    for ( ; !w.full_circle(); w = w.next().next().opp() )
    {
        Vec3d from_vertex   = m.pos( w.prev().vertex()),
              to_vertex     = m.pos( w.vertex());
        Vec3d dir           = to_vertex - from_vertex;
        points_dirs.push_back( make_pair( from_vertex, dir ));
    }
    
    for (int i = 1; i < slices; ++i)
    {
        inserted_vertices_position.push_back( vector<Vec3d>( ));
        for( point_dir pd : points_dirs )
        {
            inserted_vertices_position[i-1].push_back( pd.first + pd.second * ((double)i/(double)slices) );
        }
    }
    
    
    
    for (int i = 1; i < slices; ++i)
    {
        inserted_vertices.clear();
        Walker w = m.walker(h);
        int count = 0;
        // split all the edges and keep track of all the faces that must be split
        do
        {
            VertexID new_v = m.split_edge(w.halfedge());
            inserted_vertices.push_back( new_v );
            //update new_v position
            m.pos(new_v) = inserted_vertices_position[i-1][count];
            
            faces_to_split.push_back(w.face());
            w = w.next().next().opp();
            
            count++;
        }while ( !w.full_circle( ));
        
        // split the faces
        for (int i = 0; i < inserted_vertices.size() - 1; i++)
        {
            m.split_face_by_edge(faces_to_split[i], inserted_vertices[i], inserted_vertices[i+1]);
        }
        // split the last face ( the one that closes the loop
        m.split_face_by_edge(faces_to_split.back(), inserted_vertices[inserted_vertices.size()-1], inserted_vertices.front());

    }
    
    

    
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
        length = mean_length_of_outoing_he( m, v );
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


void scale_ring_radius ( HMesh::Manifold& m, HMesh::HalfEdgeID h, double ratio, bool smoothing )
{
    vector< VertexID > vertices;
    Vec3d barycenter = ring_barycenter( m, h, vertices );
    if( smoothing )
    {
        double max_radius = -1;
        for( auto v : vertices)
        {
            double curr_radius = ( m.pos( v ) - barycenter ).length();
            if( curr_radius > max_radius ) { max_radius = curr_radius; }
        }
        
        for( auto v : vertices)
        {
            auto dir = m.pos(v) - barycenter;
            double dir_len = dir.length();
            dir.normalize();
            m.pos(v) = barycenter + dir * (( dir_len + max_radius ) / 2 ) * ratio;
        }
    }
    else
    {
        for( auto v : vertices)
        {
            m.pos(v) = barycenter + ( m.pos( v ) - barycenter ) * ratio;
        }   
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
    vector< VertexID > neighbors;
    Walker w = m.walker(pole);
    for( ; !w.full_circle(); w = w.circulate_vertex_ccw()) { neighbors.push_back(w.vertex()); }
    // find the combination of three that creates the biggest triangle
    // take those 3 points as plane - only an approximation ( maybe for now )
    int pace = neighbors.size() / 3;
    Vec3d p0 = m.pos( neighbors[0] );
    Vec3d p1 = m.pos( neighbors[pace] );
    Vec3d p2 = ( pace * 3 < neighbors.size( ) && neighbors.size() > 4 )
               ? m.pos(neighbors[pace * 2 + pace / 2])
               : m.pos(neighbors[pace * 2]);
    // move the pole in its projection onto the plane.
    Procedural::Geometry::Plane p(p0, p1, p2);
    Vec3d new_pole = p.ortho(p0, m.pos(pole));
    m.pos(pole) = new_pole;
    
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
            
void smooth_pole ( Manifold& m, VertexID pole, HalfEdgeAttributeVector<EdgeInfo> edge_info )
{
    if ( !is_pole(m, pole) ) return;
    vector< VertexID > selected;
    Walker w = m.walker( pole );
    // found the start of the pole

//#ifdef DEBUG
//    for( HalfEdgeID he : m.halfedges()) { assert( edge_info[he].is_junction() == edge_info[ m.walker(he).opp().halfedge()].is_junction()); }
//    for( HalfEdgeID he : m.halfedges()) { assert( !edge_info[he].is_junction() ); }
//#endif
    
    // this sould be done better, smoothing all the branch, I guess
    while ( !edge_info[ w.next().halfedge() ].is_junction() )
    {
        w = w.next().opp().next();
        Walker ring = m.walker( w.next().halfedge() );
        while( !ring.full_circle() ) { selected.push_back( ring.vertex( ));
            ring = ring.next().opp().next(); }

    }
    Walker ring = m.walker( w.next().halfedge() );
    while( !ring.full_circle() ) { selected.push_back( ring.vertex( ));
                                   ring = ring.next().opp().next(); }
    
    selected_vertices_cotangent_weights_laplacian( m, selected );
}
    

}}}
