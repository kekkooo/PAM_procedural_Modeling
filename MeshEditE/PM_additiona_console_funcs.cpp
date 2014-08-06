//
//  PM_additiona_console_funcs.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 24/07/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "PM_additiona_console_funcs.h"
#include <GEL/GLGraphics/MeshEditor.h>
#include <strstream>
#include <istream>
#include <fstream>
#include <string>
#include <regex>
#include <set>
#include <queue>
#include "Test.h"
#include "polarize.h"

using namespace GLGraphics;

using namespace std;
using namespace HMesh;

// inverse distance laplacian smoothing
void console_test_inverse_distance_smoothing( MeshEditor *me, const std::vector< std::string > &args )
{
    int times=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> times;
    }
    for (int iter = 0; iter < times; iter++) {
        inverse_distance_laplacian_smoothing( me->active_mesh( ));
    }

}

// cotangent weights laplacian smoothing
void console_test_cotangent_smoothing( MeshEditor *me, const std::vector< std::string > &args )
{
    int times=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> times;
    }
    for (int iter = 0; iter < times; iter++) {
        cotangent_weights_laplacian_smoothing( me->active_mesh( ));
    }
    
}




// THE LOGIC OF THIS THING SHOULD BE MOVED INTO ANOTHER CLASS
void console_test_add_branch( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    
    Manifold& m = me->active_mesh();
    auto selection = me->get_vertex_selection();
    HMesh::VertexAttributeVector<int> ring( selection.size(), 0);
    VertexID used_vertex = InvalidVertexID;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
//    for( VertexIDIterator vit = m.vertices_begin();
//        vit != m.vertices_end() && used_vertex == InvalidVertexID; ++vit)
    for( VertexIDIterator vit = m.vertices_begin();
        vit != m.vertices_end(); ++vit)    
    {
        VertexID vid = *vit;
        map< HMesh::VertexID, CGLA::Vec3d > vert_pos;
        // need to select each vertex of the face ring of the selected vertex
        if( selection[vid] )
        {
            Walker w = m.walker(vid);
            for (; !w.full_circle(); w=w.circulate_vertex_ccw() )
            {
                if ( edge_info[w.halfedge()].is_rib())
                {
                    ring[w.vertex()] = 1;
                    ring[w.prev().vertex()] = 1;
                }
                vert_pos[ w.vertex()] = m.pos(w.vertex());
            }
            
            // if the number of outgoint halfedges >4, it was selected a pole
            if( w.no_steps() != 4 ) { return; }
            
//            vert_pos[vid] = m.pos( vid );
            
            polar_add_branch(m, ring);
            for( auto kv : vert_pos )
            {
                m.pos( kv.first ) = kv.second;
            }
            
            used_vertex = vid;
        }
    }
    
    me->active_visobj().clear_selection();
    m.cleanup();

    cout << "add branches to selected vertices" << endl;

}

void console_test_cut_branch( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    bool done = false;
    Manifold& m = me->active_mesh();
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    VertexID selected_pole = InvalidVertexID;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    HalfEdgeID selected_he = InvalidHalfEdgeID;
    
    // consider only selected vertices that are not poles
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit] )
        {
            selected.push_back(*vit);
        }
    }
    
    // take a couple of selected vertices that define a spine edge
    //    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end() && !done; ++vit)
    for( vertexID_iter vit = selected.begin();
         vit != selected.end() && ( !done || selected_pole == InvalidVertexID ); ++vit)
    {
        
        if( is_pole(m, *vit ))
        {
            selected_pole = *vit;
        }
        else if( !done )
        {
            Walker w = m.walker( *vit );
            for (; !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                assert(*vit != w.vertex());
                // no need to check if it is not a pole because from a pole you can't have ribs
                if( me->get_vertex_selection()[w.vertex()] && edge_info[w.halfedge()].is_rib()  )
                {
                    selected_he = w.halfedge();
                    done = true;
                }
            }
        }
    }
    assert(selected_he != InvalidHalfEdgeID);
    assert( selected_pole != InvalidVertexID);
//    split_ring_of_quads(m, selected_he);
    cut_branch(m, selected_he, selected_pole);
    m.cleanup();
}

// this works if you have selected 2 vertices and both of them are not poles,
void console_test_split_quad_ring( MeshEditor *me, const std::vector< std::string > &args )
{
    // for now this is useless. But it is definitely something that has to be done.
    // if number is power of 2 you can just call iteratively the same function
    // as in a mergesort.
    // if the number is not power of two, you need to operate differentely.
    // maybe you can create 2 library calls.
    int times = 1;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> times;
    }
    
    me->save_active_mesh();
    bool done = false;
    Manifold& m = me->active_mesh();
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    HalfEdgeID selected_he = InvalidHalfEdgeID;
    
    // consider only selected vertices that are not poles
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit] && !is_pole(m, *vit))
        {
            selected.push_back(*vit);
        }
    }
    
    // take a couple of selected vertices that define a spine edge
//    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end() && !done; ++vit)
    for( vertexID_iter vit = selected.begin(); vit != selected.end() && !done; ++vit)
    {
        Walker w = m.walker( *vit );
        for (; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            assert(*vit != w.vertex());
            if( me->get_vertex_selection()[w.vertex()] && edge_info[w.halfedge()].is_spine()  )
            {
                selected_he = w.halfedge();
                done = true;
            }
        }
    }
    split_ring_of_quads(m, selected_he);
    me->active_visobj().clear_selection();
    m.cleanup();
}

// this works if you have selected 2 vertices
void console_test_split_pole_to_pole( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    bool done = false;
    Manifold& m = me->active_mesh();
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
    HalfEdgeID selected_he = InvalidHalfEdgeID;
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            selected.push_back(*vit);
        }
    }
    
    // take a couple of selected vertices that define a rib edge
    for( vertexID_iter vit = selected.begin(); vit != selected.end() && !done; ++vit)
    {
        Walker w = m.walker( *vit );
        for (; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            assert(*vit != w.vertex());
            if( me->get_vertex_selection()[w.vertex()] && edge_info[w.halfedge()].is_rib()  )
            {
                selected_he = w.halfedge();
                done = true;
            }
        }
    }
    
    split_from_pole_to_pole( m, selected_he );
    me->active_visobj().clear_selection();
    m.cleanup();
}

void console_test_extrude_vertices( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    double length = NAN;
    bool add_quads_on_poles = false;
    double ratio = 1.0;

    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> length;
    }
    
    if(args.size() > 1){
        istringstream a1(args[1]);
        a1 >> add_quads_on_poles;
    }
    
    if(args.size() > 2){
        istringstream a2(args[2]);
        a2 >> ratio;
    }
    
    cout << length << " # " << add_quads_on_poles << endl;
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            extrude_pole(m, *vit, length, add_quads_on_poles, ratio);
        }
    }
    m.cleanup();
}


void console_test_extrude_vertices_alt( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    double length = NAN;
    bool add_quads_on_poles = false;
    double ratio = 1.0;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> length;
    }
    
    if(args.size() > 1){
        istringstream a1(args[1]);
        a1 >> add_quads_on_poles;
    }
    
    if(args.size() > 2){
        istringstream a2(args[2]);
        a2 >> ratio;
    }
    
    cout << length << " # " << add_quads_on_poles;// << endl;
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {

            if( is_pole( m, *vit ))
            {
                // select direction
                Walker w = m.walker(*vit);
                for( ; !w.full_circle() && me->get_vertex_selection()[w.vertex()] != 1; w = w.circulate_vertex_ccw() );

                CGLA::Vec3d dir = m.pos(*vit) - m.pos(w.vertex());
                dir.normalize();
                dir = dir * length;

                extrude_pole(m, *vit, dir * length, add_quads_on_poles, ratio);
            }
        }
    }
    m.cleanup();
}




void console_test_set_ring_radius( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();

    double radius = 1.0;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> radius;
    }
    
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
    HalfEdgeID selected_he = InvalidHalfEdgeID;
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end()  && selected.size() < 2; ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            selected.push_back(*vit);
        }
    }
    
    selected_he = find_rib_edge(m, selected[0], selected[1], edge_info);
    if( selected_he != InvalidHalfEdgeID)
    {
        set_ring_radius(m, selected_he, radius);
    }
    m.cleanup();
}

void console_test_scale_ring_radius( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    
    double ratio = 1.0;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> ratio;
    }
    
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
    HalfEdgeID selected_he = InvalidHalfEdgeID;
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end()  && selected.size() < 2; ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            selected.push_back(*vit);
        }
    }
    
    selected_he = find_rib_edge(m, selected[0], selected[1], edge_info);
    if( selected_he != InvalidHalfEdgeID)
    {
//        set_ring_radius(m, selected_he, radius);
        scale_ring_radius(m, selected_he, ratio);
    }
    m.cleanup();
}

void console_test_scale_selected_rings( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    
    vector< HalfEdgeID > selected_rings;
    
    double ratio = 1.0;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> ratio;
    }
    
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            selected.push_back(*vit);
        }
    }
    
    for( VertexID vid : selected )
    {
        Walker w = m.walker(vid);
        for( ; !w.full_circle() && !edge_info[w.halfedge()].is_rib(); w = w.circulate_vertex_ccw() );
        if( edge_info[w.halfedge()].is_rib() ) selected_rings.push_back(w.halfedge());
    }
    

    for( auto selected_he : selected_rings )
    {
        scale_ring_radius(m, selected_he, ratio);
    }
    m.cleanup();
}

void console_test_perturbate( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    
    double cutoff   = 0.1;
    double ratio    = 1.0;
    int type        = 0;
    int vtype       = VertexType::POLE;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> type;
    }
    
    if(args.size() > 1){
        istringstream a1(args[1]);
        a1 >> ratio;
    }
    
    if(args.size() > 2){
        istringstream a1(args[2]);
        a1 >> cutoff;
    }

    
    if (type == 1) { vtype = VertexType::REGULAR; }
    if( type == 2) { vtype |= VertexType::REGULAR; }

    cout << " perturbating " << vtype << " with ratio : " << ratio;
    add_noise(m, vtype, ratio, cutoff);
    
}



void register_more_basic_console_funcs( GLGraphics::MeshEditor *me )
{
    me->register_console_function( "test.smoothing.inverse_distance",
                                  console_test_inverse_distance_smoothing,
                                   "test.smoothing.inverse_distance <iterations>" );
    
    me->register_console_function( "test.smoothing.cotangent",
                                  console_test_cotangent_smoothing,
                                  "test.smoothing.cotangent <iterations>" );


    
    me->register_console_function( "test.add_branch", console_test_add_branch,
                                  "test.add_branch" );
    
    me->register_console_function( "test.cut_branch", console_test_cut_branch,
                                  "test.cut_branch" );

    me->register_console_function( "test.split_quad_ring", console_test_split_quad_ring,
                                  "test.split_quad_ring" );

    me->register_console_function( "test.split_pole_to_pole", console_test_split_pole_to_pole,
                                  "test.split_pole_to_pole" );
    
    me->register_console_function( "test.extrude_poles", console_test_extrude_vertices,
                                  "test.extrude_poles" );

    me->register_console_function( "test.extrude_poles_alt", console_test_extrude_vertices_alt,
                                  "test.extrude_poles_alt" );
    
    me->register_console_function( "test.set_ring_radius", console_test_set_ring_radius,
                                  "test.set_ring_radius" );
    
    me->register_console_function( "test.scale_ring_radius", console_test_scale_ring_radius,
                                  "test.scale_ring_radius" );
    
    me->register_console_function( "test.perturbate", console_test_perturbate,
                                  "test.perturbate" );
    
    me->register_console_function( "test.scale_selected_rings", console_test_scale_selected_rings,
                                  "test.scale_selected_rings" );
}