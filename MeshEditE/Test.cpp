//
//  Test.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 25/07/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "Test.h"
#include <GEL/GLGraphics/MeshEditor.h>
#include <GEL/CGLA/Vec3d.h>

#include <map>
#include <queue>
#include <algorithm>

#include "Procedural/Helpers/Plane.h"

using namespace GLGraphics;
using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace Procedural::Geometry;

void create_cube_vertices( double edge_length, vector<Vec3f> &vertices )
{
    float plus  = edge_length / 2.0;
    float minus = -plus;

    // initialize points
    vertices[0] = Vec3f( minus, plus, plus );
    vertices[1] = Vec3f( minus, plus, minus );
    vertices[2] = Vec3f( plus, plus, minus );
    vertices[3] = Vec3f( plus, plus, plus );
    vertices[4] = Vec3f( minus, minus, plus );
    vertices[5] = Vec3f( minus, minus, minus );
    vertices[6] = Vec3f( plus, minus, minus );
    vertices[7] = Vec3f( plus, minus, plus );
}

// Creates a cube with little piramids (with squared foot) instead of top and bottom face
// One of the simplest PAMs
void create_basic_PAM( HMesh::Manifold& m, double edge_length )
{
    vector< Vec3f > vertices( 10 );
    vector< int >   faces   ( 12, 3 );
    vector< int >   indices ( 40 );
    
    create_cube_vertices(edge_length, vertices);
    // add the pole vertices;
    vertices[8] = Vec3f(0.0, ( 2.0 * edge_length / 3.0 ), 0.0 );
    vertices[9] = Vec3f(0.0, ( -2.0 * edge_length / 3.0 ), 0.0 );
    
    faces[0] = 4;   faces[1] = 4;   faces[2] = 4;   faces[3] = 4;
    
    // LEFT
    indices[0]  = 1;    indices[1]  = 5;    indices[2]  = 4;    indices[3]  = 0;
    // BACK
    indices[4]  = 1;    indices[5]  = 2;    indices[6]  = 6;    indices[7]  = 5;
    // FRONT
    indices[8]  = 3;    indices[9]  = 0;    indices[10] = 4;    indices[11] = 7;
    // RIGHT
    indices[12] = 7;    indices[13] = 6;    indices[14] = 2;    indices[15] = 3;
    // upper pole
    indices[16] = 2;    indices[17] = 1;    indices[18] = 8;
    indices[19] = 1;    indices[20] = 0;    indices[21] = 8;
    indices[22] = 0;    indices[23] = 3;    indices[24] = 8;
    indices[25] = 3;    indices[26] = 2;    indices[27] = 8;
    // lower pole
    indices[28] = 7;    indices[29] = 4;    indices[30] = 9;
    indices[31] = 4;    indices[32] = 5;    indices[33] = 9;
    indices[34] = 5;    indices[35] = 6;    indices[36] = 9;
    indices[37] = 6;    indices[38] = 7;    indices[39] = 9;
    
    m.build( vertices.size( ), reinterpret_cast< float * >( &vertices[0] ),
            faces.size( ), &faces[0], &indices[0] );
}

// Creates the triangle mesh of a cube
void create_triangle_cube( HMesh::Manifold& m, double edge_length )
{
    vector< Vec3f > vertices( 8 );
    vector< int >   faces   ( 12, 3 );
    vector< int >   indices ( 36 );
    
    create_cube_vertices(edge_length, vertices);
    
    // for each face save the indices of its vertices
    // UP
    // face 0
    indices[0] = 2;     indices[1] = 1;     indices[2] = 3;
    // face 1
    indices[3] = 3;     indices[4] = 1;     indices[5] = 0;
    // LEFT
    // face 2
    indices[6] = 1;     indices[7] = 5;     indices[8] = 0;
    // face 3
    indices[9] = 0;     indices[10] = 5;    indices[11] = 4;
    // BACK
    // face 4
    indices[12] = 1;    indices[13] = 2;    indices[14] = 5;
    // face 5
    indices[15] = 5;    indices[16] = 2;    indices[17] = 6;
    // FRONT
    // face 6
    indices[18] = 3;    indices[19] = 0;    indices[20] = 7;
    // face 7
    indices[21] = 7;    indices[22] = 0;    indices[23] = 4;
    // RIGHT
    // RIGHT 8
    indices[24] = 7;    indices[25] = 2;    indices[26] = 3;
    // face 9
    indices[27] = 2;    indices[28] = 7;    indices[29] = 6;
    // BOTTOM
    // face 10
    indices[30] = 7;    indices[31] = 5;    indices[32] = 6;
    // face 11
    indices[33] = 5;    indices[34] = 7;    indices[35] = 4;
    
    m.build( vertices.size( ), reinterpret_cast< float * >( &vertices[0] ),
            faces.size( ), &faces[0], &indices[0] );
}

// Creates the quad mesh of a cube
void create_quads_cube( HMesh::Manifold& m, double edge_length )
{
    vector< Vec3f > vertices( 8 );
    vector< int >   faces   ( 6, 4 );
    vector< int >   indices ( 24 );
    
    create_cube_vertices(edge_length, vertices);
   
    // for each face save the indices of its vertices
    // UP
    indices[0]  = 2;    indices[1]  = 1;    indices[2]  = 0;    indices[3]  = 3;
    // LEFT
    indices[4]  = 1;    indices[5]  = 5;    indices[6]  = 4;    indices[7]  = 0;
    // BACK
    indices[8]  = 1;    indices[9]  = 2;    indices[10] = 6;    indices[11] = 5;
    // FRONT
    indices[12] = 3;    indices[13] = 0;    indices[14] = 4;    indices[15] = 7;
    // RIGHT
    indices[16] = 7;    indices[17] = 6;    indices[18] = 2;    indices[19] = 3;
    // BOTTOM
    indices[20] = 7;    indices[21] = 4;    indices[22] = 5;    indices[23] = 6;
    
    m.build( vertices.size( ), reinterpret_cast< float * >( &vertices[0] ),
            faces.size( ), &faces[0], &indices[0] );
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
        add_ring_around_pole(m, v, scaling);
    }
}

void add_ring_around_pole ( Manifold& m, VertexID pole, double scaling )
{
    typedef vector< pair< VertexID, VertexID > >        NeighborVector;
    
    auto new_vs = split_ring_of_quads( m, m.walker( pole ).halfedge( ));
    // if it is selected to preserve the radius of the shape in the inserction point of the triangles fan
    // THIS COULD BE DONE WITH the function to scale the ring.
    if( scaling > 0 )
    {
        NeighborVector nb;
        HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
        // build corresponnces between the newly inserted vertices and
        // the next in the spine edge's vertex chain (nsevc vertex)
        for( VertexID new_v : new_vs )
        {
            Walker w = m.walker(new_v);
            bool done = false;
            // must select the correct half edge, not the one that has the pole as its endpoint
            for (; !w.full_circle() && !done; w = w.circulate_vertex_ccw())
            {
                if ( edge_info[w.halfedge()].is_spine() && w.vertex() != pole )
                {
                    nb.push_back( make_pair( new_v, w.vertex()));
                    done = true;
                }
            }
        }
        // find barycenter and distances from it for the nsevc vertices ring
        Vec3d nsevc_barycenter(0), barycenter(0);
        for( auto corr : nb) {  barycenter += m.pos(corr.first);
                                nsevc_barycenter += m.pos(corr.second); }
        barycenter       /= nb.size();
        nsevc_barycenter /= nb.size();
        // this vector will store the distance of the nsevc vertex from the barycenter of its ring,
        // but bonding it to the correspondent new vertex, also storing the direction in which the vertices must move
        DistanceDirectionVector dists;
        for( auto corr : nb)
        {
            assert( corr.first != pole && corr.second != pole );

            Vec3d dir = ( m.pos( corr.first ) - barycenter );
            dir.normalize();
            dists.push_back( make_tuple( corr.first, (m.pos(corr.second) - nsevc_barycenter).length( ), dir ));
        }
        
        // update position according to original direction, previous ring distances from barycenter and preserving ratio
        for( auto triple : dists )
        {
         //   cout << " dist " << get<1>( triple ) << endl;
            m.pos( get<0>( triple )) = barycenter + ( get<2>( triple ) * get<1>( triple ) * scaling);
        }
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
    //only for debug
    test_ring_barycenter(m, h);
    // end - Only for debug
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

void cut_branch( Manifold& m, HalfEdgeID h, VertexID pole )
{
//  0) find the ring of vertices
    vector< VertexID > vertices;

    // 1) save a vector from one of the pole's neighbors to the pole ( or I need all of them? )
    // need to save anche which of the vertices of the slitted ring has to be used as reference
    Vec3d       barycenter      = ring_barycenter( m, h, vertices );
    Vec3d       pole_dir        = m.pos(pole) - barycenter;
    // 2) slit that ring of vertices
    VertexAttributeVector<int> selection(m.no_vertices(), 0);
    for( VertexID v : vertices ) { selection[v] = 1; }
    HalfEdgeID hole_edge = m.slit_edges(selection);
    // 3) delete all the faces of the slitted component
    // this should cut also the sub_branches!!!!!!!!!!!! ( NOT SO EASY! )
    vector< VertexID > v_to_delete;
    v_to_delete.push_back(pole);
    Walker  w       = m.walker(pole);
    bool    done    = false;
            w       = m.walker( w.next().halfedge( ));
    
    while( !done )
    {
        // this just works for a single branch with no sub-branches
        /* starting from the pole
            for each outgoint branch
                mark the vertex to be deleted and move the walker to next.opp.next
            // each walker stops when reaches a vertex inside "vertices" ( put them in a SET! )         
         */
        
        
        
        
//        done = ( w.opp().face() == InvalidFaceID);
//        for (; !w.full_circle(); w = w.next().opp().next()) {
//            v_to_delete.push_back(w.vertex());
//        }
//        if (!done)
//        {
//            w = m.walker( w.opp().next().next().halfedge( ));
//        }
    }
    for ( auto vtd : v_to_delete )
    {
        cout << "deleting vertex " << vtd << endl;
        m.remove_vertex(vtd);
    }
    // 4) rebuild the pole
    FaceID      f           = m.close_hole( hole_edge );
    VertexID    new_pole    = m.split_face_by_vertex( f );
    Vec3d       f_centroid  = ring_barycenter( m, hole_edge );
    m.pos( new_pole )       = f_centroid + pole_dir;
}

// this works only for cutting branches that are not part of a loop
// NEED TO KNOW how to understand if the ring is part of a loop
// ALSO need to understend in which direction should the walker move in order to reach the pole
void cut_branch ( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    // if not loop OR not between two joints
    //              find the right pole
    //              call cut_branch( m, h, pole )
    // else cut creating two poles with a simple policy
}



/* THESE ARE UTILITIES AND SHOULD BE MOVED SOMEWHERE ELSE */
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


bool test_ring_barycenter( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    vector< VertexID > vs;
    Vec3d b = ring_barycenter(m, h, vs);
    assert( vs.size() > 0 );
    Plane p( m.pos( vs.front()), m.pos(vs.back()), m.pos(vs[ vs.size()/2 ]));
    
    for( VertexID v : vs )
    {
        assert( p.OnPlane( m.pos( v )));
    }
    assert( p.OnPlane( b));

    return true;
}

CGLA::Vec3d ring_barycenter ( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    vector< VertexID > vs;
    return ring_barycenter(m, h, vs);
}

CGLA::Vec3d ring_barycenter ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                              std::vector< HMesh::VertexID > &vertices   )
{
    assert(h != InvalidHalfEdgeID);
    
    Walker w = m.walker( h );
    Vec3d barycenter(0);
    
    while ( !w.full_circle( ))
    {
        vertices.push_back( w.vertex( ));
        barycenter += m.pos( w.vertex() );
        w = w.next().opp().next();
    }
    
    barycenter /= vertices.size();
    return barycenter;
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


/* THESE ARE UTILITIES AND SHOULD BE MOVED SOMEWHERE ELSE */

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

Vec3d simple_random_direction ( Manifold& m, VertexID v, double vertex_normal_weight )
{
    // take the vertex normal
    // choose randomly the normal of one of the sorrounding faces
    // sum them
    // multiply for a random ration between the range [ -0.5, 0.5 ]
    Vec3d v_normal = vertex_normal( m, v );
    int     times   = rand() % 23;
    double  ratio   = ( (double)( rand() % 1000 ) - 500.0 )/10000.0;
    cout << "times is : " << times << " ratio is : " << ratio <<endl;
    Walker w = m.walker(v);
    for( int i = 0; i < times; w = w.circulate_vertex_ccw(), i++ );
    
    Vec3d f_normal  = face_normal( m, w.face( ));
    
    f_normal.normalize();
    v_normal.normalize();
    
    Vec3d dir       = ( v_normal + f_normal ) * ratio;
    return dir;
}

Vec3d alt_simple_random_direction ( Manifold& m, VertexID   v )
{
    Vec3d v_normal = vertex_normal( m, v );
    int     times   = rand() % 23;
    double  ratio   = ( (double)( rand() % 1000 ) - 500.0 )/10000.0;
    cout << "times is : " << times << " ratio is : " << ratio <<endl;
    Walker w = m.walker(v);
    for( int i = 0; i < times; w = w.circulate_vertex_ccw(), i++ );
    
    Vec3d other_dir  = m.pos( w.vertex()) - m.pos(v);
    
    other_dir.normalize();
    v_normal.normalize();
    
    Vec3d dir       = ( v_normal + other_dir * 1.5 ) * ratio;
    return dir;

}

