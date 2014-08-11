//
//  geometric_console_functions.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "geometric_console_functions.h"

#include <GEL/GLGraphics/MeshEditor.h>
#include <MeshEditE/Procedural/Operations/geometric_operations.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <strstream>
#include <istream>
#include <fstream>
#include <string>
#include <regex>
#include <set>
#include <queue>
#include "polarize.h"

using namespace GLGraphics;
using namespace std;
using namespace HMesh;
using namespace Procedural::Operations::Geometric;
using namespace Procedural::Structure;

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
    bool        done        = false;
    Manifold&   m           = me->active_mesh();
    HalfEdgeID  selected_he = InvalidHalfEdgeID;
    vector< VertexID > selected;
    typedef vector< VertexID >::iterator vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );

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
    for( vertexID_iter vit = selected.begin(); vit != selected.end() && !done; ++vit )
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
    split_ring_of_quads( m, selected_he );
    me->active_visobj().clear_selection();
    m.cleanup();
}

// this works if you have selected 2 vertices
void console_test_split_pole_to_pole( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    bool        done        = false;
    Manifold&   m           = me->active_mesh();
    HalfEdgeID selected_he  = InvalidHalfEdgeID;

    HalfEdgeAttributeVector<EdgeInfo>
                edge_info   = label_PAM_edges( m );

    vector< VertexID >                      selected;
    typedef vector< VertexID >::iterator    vertexID_iter;
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            selected.push_back(*vit);
        }
    }
    
    // take a couple of selected vertices that define a rib edge
    for( vertexID_iter vit = selected.begin(); vit != selected.end() && !done; ++vit )
    {
        Walker w = m.walker( *vit );
        for (; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            assert(*vit != w.vertex());
            if( me->get_vertex_selection()[w.vertex()] && edge_info[w.halfedge()].is_rib()  )
            {
                selected_he = w.halfedge();
                done        = true;
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
    Manifold&   m                   = me->active_mesh();
    double      length              = NAN;
    bool        add_quads_on_poles  = false;
    double      ratio               = 1.0;
    
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
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if ( me->get_vertex_selection()[*vit] )
        {
            extrude_pole( m, *vit, length, add_quads_on_poles, ratio );
        }
    }
    m.cleanup();
}

// QUESTA ROBA VA SISTEMATO DECENTEMENTE.
void console_test_extrude_vertices_alt( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m                   = me->active_mesh();
    double      length              = NAN;
    bool        add_quads_on_poles  = false;
    double      ratio               = 1.0;
    
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
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit )
    {
        if( me->get_vertex_selection()[*vit] )
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
    Manifold&   m           = me->active_mesh();
    double      radius      = 1.0;
    HalfEdgeID  selected_he = InvalidHalfEdgeID;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> radius;
    }
    
    vector< VertexID >                      selected;
    typedef vector< VertexID >::iterator    vertexID_iter;
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
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
        set_ring_radius( m, selected_he, radius );
    }
    m.cleanup();
}

void console_test_scale_ring_radius( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m           = me->active_mesh();
    double      ratio       = 1.0;
    HalfEdgeID  selected_he = InvalidHalfEdgeID;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> ratio;
    }
    
    vector< VertexID >                      selected;
    typedef vector< VertexID >::iterator    vertexID_iter;
    HalfEdgeAttributeVector<EdgeInfo>       edge_info = label_PAM_edges( m );
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end()  && selected.size() < 2; ++vit)
    {
        if (me->get_vertex_selection()[*vit])
        {
            selected.push_back(*vit);
        }
    }
    
    selected_he = find_rib_edge( m, selected[0], selected[1], edge_info );
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
    Manifold&   m       = me->active_mesh();
    double      ratio   = 1.0;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> ratio;
    }
    
    set< HalfEdgeID >                       selected_rings;
    vector< VertexID >                      selected;
    typedef vector< VertexID >::iterator    vertexID_iter;
    HalfEdgeAttributeVector<EdgeInfo>       edge_info = label_PAM_edges( m );
    
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
        if( edge_info[w.halfedge()].is_rib() && selected_rings.count(w.opp().halfedge()) == 0 )
        {
            selected_rings.insert(w.halfedge());
        }
    }
    
    for( auto selected_he : selected_rings )
    {
        scale_ring_radius( m, selected_he, ratio );
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
    
    
    if ( type == 1 ) { vtype = VertexType::REGULAR; }
    if ( type == 2 ) { vtype |= VertexType::REGULAR; }
    
    cout << " perturbating " << vtype << " with ratio : " << ratio;
    add_noise(m, vtype, ratio, cutoff);
}

namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_geometric_console_funcs(GLGraphics::MeshEditor* me)
        {
            me->register_console_function( "test.geometry.split_quad_ring", console_test_split_quad_ring,
                                           "test.geometry.split_quad_ring" );
            
            me->register_console_function( "test.geometry.split_pole_to_pole", console_test_split_pole_to_pole,
                                           "test.geometry.split_pole_to_pole" );
            
            me->register_console_function( "test.geometry.extrude_poles", console_test_extrude_vertices,
                                           "test.extrude_poles" );
            
            me->register_console_function( "test.geometry.extrude_poles_alt", console_test_extrude_vertices_alt,
                                           "test.geometry.extrude_poles_alt" );
            
            me->register_console_function( "test.geometry.set_ring_radius", console_test_set_ring_radius,
                                           "test.geometry.set_ring_radius" );
            
            me->register_console_function( "test.geometry.scale_ring_radius", console_test_scale_ring_radius,
                                           "test.geometry.scale_ring_radius" );
            
            me->register_console_function( "test.geometry.perturbate", console_test_perturbate,
                                           "test.geometry.perturbate" );
            
            me->register_console_function( "test.geometry.scale_selected_rings", console_test_scale_selected_rings,
                                           "test.geometry.scale_selected_rings" );
        }
}}