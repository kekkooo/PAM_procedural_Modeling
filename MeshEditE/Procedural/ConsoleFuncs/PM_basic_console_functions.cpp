//
//  PM_basic_console_functions.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 06/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "PM_basic_console_functions.h"
#include <GEL/GLGraphics/MeshEditor.h>
#include <sstream>
//#include <strstream>
//#include <istream>
//#include <fstream>
//#include <string>
//#include <regex>
//#include <set>
//#include <queue>
//#include "Test.h"
#include "polarize.h"
#include <MeshEditE/Procedural/Operations/basic_shapes.h>
#include <MeshEditE/Procedural/Operations/Algorithms.h>
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>


using namespace GLGraphics;
using namespace Procedural::Operations;
using namespace Procedural::Operations::Algorithms;

using namespace std;
using namespace HMesh;

/*****************************************************************************************
**                           register_basic_console_funcs                               **
*****************************************************************************************/

// Console call to create a cube built with triangular faces
void console_test_cube_triangles( MeshEditor *me, const std::vector< std::string > &args )
{
    double side = 0.6;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> side;
    }
    cout << "triangulated cube" << endl;
    create_triangle_cube(me->active_mesh(), side);
}

// Console call to create a cube built with quadrangular faces
void console_test_cube_quads( MeshEditor *me, const std::vector< std::string > &args )
{
    double side = 0.6;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> side;
    }
    cout << "quadrangulated cube" << endl;
    create_quads_cube(me->active_mesh(), side);
}

// Console call to create a cube with poles replacing top and bottom faces
void console_test_basic_PAM( MeshEditor *me, const std::vector< std::string > &args )
{
    double side = 0.6;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> side;
    }
    cout << "basic pam" << endl;
    create_basic_PAM(me->active_mesh(), side);
}

/*****************************************************************************************
**                         register_algorithm_console_funcs                             **
*****************************************************************************************/

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
    me->save_active_mesh();
    Manifold&   m       = me->active_mesh();
    int         times   = 1;
    vector< VertexID > selected;

    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> times;
    }
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit )
    {
        if( me->get_vertex_selection()[*vit] ) selected.push_back(*vit);
    }
    for (int iter = 0; iter < times; iter++) {
//        cotangent_weights_laplacian_smoothing( me->active_mesh( ));
        selected_vertices_cotangent_weights_laplacian( m, selected );
    }
}

// inverse distance laplacian smoothing
void console_test_selected_inverse_distance_smoothing( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m       = me->active_mesh();
    int         times   = 1;
    vector< VertexID > selected;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> times;
    }
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit )
    {
        if( me->get_vertex_selection()[*vit] ) selected.push_back(*vit);
    }
    
    for (int iter = 0; iter < times; iter++) {
        selected_vertices_inverse_distance_laplacian( m, selected );
//        inverse_distance_laplacian_smoothing( me->active_mesh( ));
    }
    

}

// cotangent weights laplacian smoothing
void console_test_selected_cotangent_smoothing( MeshEditor *me, const std::vector< std::string > &args )
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


void console_test_slit_selection( MeshEditor *me, const std::vector< std::string > &args )
{
    me->active_mesh().slit_edges(me->get_vertex_selection());
}


void console_test_delete_vertices( MeshEditor *me, const std::vector< std::string > &args )
{
    auto *m = &(me->active_mesh());
    for( VertexIDIterator vit = m->vertices_begin(); vit != m->vertices_end(); ++vit )
    {
        if( me->get_vertex_selection()[*vit] )
            m->remove_vertex(*vit);
    }
    m->cleanup();
}

void console_test_print_vertex_info( MeshEditor *me, const std::vector< std::string > &args )
{
    auto &m = me->active_mesh();
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    VertexAttributeVector<Procedural::Geometry::DistanceMetrics> pole_dist;
    VertexAttributeVector<Procedural::Geometry::DistanceMetrics> junction_dist;
    VertexAttributeVector<Procedural::Geometry::DistanceMetrics> combined_dist;
    VertexAttributeVector<double> angles;
    VertexAttributeVector<double> spine_angles;
    VertexAttributeVector<double> rib_angles;
    
    Procedural::Structure::LabelJunctions( m, edge_info );
    
    Procedural::Geometry::dihedral_angles( m, edge_info, angles );
    Procedural::Geometry::dihedral_angles( m, edge_info, spine_angles, SPINE );
    Procedural::Geometry::dihedral_angles( m, edge_info, rib_angles, RIB );

    Procedural::Geometry::distance_from_poles( m, edge_info, pole_dist );
    Procedural::Geometry::distance_from_junctions( m, edge_info, junction_dist );
    Procedural::Geometry::distance_from_poles_and_junctions( m, edge_info, combined_dist );

    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit )
    {
        std::ostringstream oss;
        if( me->get_vertex_selection()[*vit] == 1 )
        {
            oss << " vertex ID                  : " << *vit                        << endl;
            oss << " is pole?                   : " << is_pole( m, *vit )          << endl;
            oss << " valence                    : " << valency( m, *vit )          << endl;
            if( !is_pole(m, *vit))
            {
            oss << " distance from pole         : " << pole_dist[*vit].first       << endl;
            oss << " distance from junction     : " << junction_dist[*vit].first   << endl;
            oss << " combined distance          : " << combined_dist[*vit].first   << endl;
            oss << " spines dihedral angle      : " << spine_angles[*vit]          << endl;
            oss << " ribs dihedral angle        : " << rib_angles[*vit]            << endl;
            oss << " combined dihedral angle    : " << angles[*vit]                << endl;
            }
            me->printf(oss.str().c_str());      oss.clear();
        }
    }
}


namespace Procedural{
    namespace ConsoleFuncs{

void register_basic_console_funcs(GLGraphics::MeshEditor* me)
{
    me->register_console_function( "test.shapes.cube_t",     console_test_cube_triangles,    "test.shapes.cube_t"     );
    me->register_console_function( "test.shapes.cube_q",     console_test_cube_quads,        "test.shapes.cube_q"     );
    me->register_console_function( "test.shapes.basic_PAM",  console_test_basic_PAM,         "test.shapes.basic_PAM"  );

}

void register_algorithm_console_funcs(GLGraphics::MeshEditor* me)
{
    me->register_console_function(
        "test.smoothing.inverse_distance",
        console_test_inverse_distance_smoothing,
        "test.smoothing.inverse_distance <iterations>"          );
    
    me->register_console_function(
        "test.smoothing.cotangent",
        console_test_cotangent_smoothing,
        "test.smoothing.cotangent <iterations>"                 );

    me->register_console_function(
        "test.smoothing.selected_inverse_distance",
        console_test_inverse_distance_smoothing,
        "test.smoothing.selected_inverse_distance <iterations>" );
    
    me->register_console_function(
        "test.smoothing.selected_cotangent",
        console_test_cotangent_smoothing,
        "test.smoothing.selected_cotangent <iterations>"        );
        
    me->register_console_function(
      "test.slit_selection",
      console_test_slit_selection,
      "test.slit_selection"        );
    
    me->register_console_function(
      "test.delete_vertices",
      console_test_delete_vertices,
      "test.delete_vertices"        );
    me->register_console_function(
      "test.print_vertex_info",
      console_test_print_vertex_info,
      "test.print_vertex_info"        );

    
    



}

}}