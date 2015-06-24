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
#include <random>
#include <GEL/HMesh/obj_save.h>
#include "polarize.h"
#include <MeshEditE/Procedural/Operations/basic_shapes.h>
#include <MeshEditE/Procedural/Operations/Algorithms.h>
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/PMEngine.h>
#include <MeshEditE/Procedural/Matches/Matches.h>
#include "patch_mapping.h"
#include <MeshEditE/Test.h>
#include <MeshEditE/Procedural/Helpers/manifold_copy.h>
#include <MeshEditE/Tests/test_random_rotation.h>


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
            oss << " distance from pole         : " << pole_dist[*vit]       << endl;
            oss << " distance from junction     : " << junction_dist[*vit]   << endl;
            oss << " combined distance          : " << combined_dist[*vit]   << endl;
            oss << " spines dihedral angle      : " << spine_angles[*vit]          << endl;
            oss << " ribs dihedral angle        : " << rib_angles[*vit]            << endl;
            oss << " combined dihedral angle    : " << angles[*vit]                << endl;
            }
            me->printf(oss.str().c_str());      oss.clear();
        }
    }
}


void console_test_save_to_results_folder( MeshEditor *me, const std::vector< std::string > &args )
{
    const string& file_name = args[0];
    if(args.size() == 1){
        if(file_name.substr(file_name.length()-4,file_name.length())==".obj")
        {
            stringstream oss;
            oss << "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/" << file_name;
            obj_save(oss.str(), me->active_mesh());
        }
    }
}


void test_match( MeshEditor *me, const std::vector< std::string > &args )
{
    static std::mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31,
                0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253> rrrr;

    rrrr.seed( time(0) );
    cout << rrrr();
    CGLA::gel_srand( rrrr() );
//    Procedural::Matching::build(me->active_mesh());
    Procedural::Matching::console_call(me->active_mesh());
}

void test_trasform( MeshEditor *me, const std::vector< std::string > &args )
{
    static std::mersenne_twister_engine<std::uint_fast32_t, 32, 624, 397, 31,
    0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253> rrrr;
    rrrr.seed( time(0) );
    cout << rrrr();
    CGLA::gel_srand( rrrr() );
    
    int id = 0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> id;
    }

    Procedural::Matching::console_call2( me->active_mesh( ), id);
}

void test_align_selected( MeshEditor *me, const std::vector< std::string > &args )
{
    Manifold& m = me->active_mesh();
    vector< VertexID > selected;
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit )
    {
        if( me->get_vertex_selection()[*vit] ) selected.push_back(*vit);
    }
    
    Procedural::Matching::align_to_selected( m, selected );
}

void console_save_patched_obj( MeshEditor *me, const std::vector< std::string > &args )
{
    Manifold& m = me->active_mesh();
    const string& file_name = args[0];
    string path;
    if(args.size() == 1){
    stringstream oss;
    oss << "/Users/francescousai/Documents/Dottorato/Visiting/" << file_name;
    path = oss.str(), me->active_mesh();

    }

    save_colored_obj( m, path );
}

void console_build_patches( MeshEditor *me, const std::vector< std::string > &args )
{
    Manifold& m = me->active_mesh();
    HMesh::FaceAttributeVector<int> face_to_patch;
    build_patches( m, face_to_patch );
}

void copy_and_delete_test( MeshEditor *me, const std::vector< std::string > &args )
{
    std::set<VertexID> m0v, m1v;
    assert( me->get_mesh(0).no_vertices() > 0 );
    assert( me->get_mesh(1).no_vertices() > 0 );
    Procedural::Helpers::add_manifold(me->get_mesh(0), me->get_mesh(1), m0v, m1v);
    Procedural::Helpers::test_delete(me->get_mesh(0), m1v );
}

void poles_info( MeshEditor *me, const std::vector< std::string > &args )
{
    for( VertexID id : me->active_mesh().vertices())
    {
        if( is_pole( me->active_mesh(), id ))
        {
            int v = valency( me->active_mesh(), id );
//            cout << id << ") has valence " << v << endl;
            me->printf("%d) has valence : %d", id, v );
        }
    }
}


/*=========================================================================*
 *                     Test Console Funcs                                  *
 *=========================================================================*/

void console_test_set_host_and_module( MeshEditor *me, const std::vector< std::string > &args )
{
    Manifold module;
    stringstream oss;
    string filename = "1.obj";
    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> filename;
    }

    oss << "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/Modules/" << filename;
    obj_load( oss.str(), module );

    Tests::RandomRotationTest &test = Tests::RandomRotationTest::getInstance();
    test.setHost( me->active_mesh() );
    test.setModule( module );
}

void console_test_randomly_rotate( MeshEditor *me, const std::vector< std::string > &args )
{
    Tests::RandomRotationTest &test = Tests::RandomRotationTest::getInstance();
    test.resetOriginalModulePosition();
    test.rotate();
}

void console_test_translate_outside_bsphere( MeshEditor *me, const std::vector< std::string > &args )
{
    Tests::RandomRotationTest &test = Tests::RandomRotationTest::getInstance();
    test.translateOutsideBSphere();
}



namespace Procedural{
    namespace ConsoleFuncs{
        
void register_match_console_funcs( GLGraphics::MeshEditor* me )
{
    me->register_console_function( "test.match.do",         test_match, "test.match.do" );
    me->register_console_function( "test.match.transform",  test_trasform, "test.match.transform" );
    me->register_console_function( "test.match.align_selected",  test_align_selected, "test.match.align_selected" );
    me->register_console_function( "test.copy_and_delete_test",  copy_and_delete_test, "copy_and_delete_test" );
    me->register_console_function( "test.poles_info",  poles_info, "test.poles_info" );
}

void register_basic_console_funcs( GLGraphics::MeshEditor* me )
{
    me->register_console_function( "test.shapes.cube_t",     console_test_cube_triangles,    "test.shapes.cube_t"     );
    me->register_console_function( "test.shapes.cube_q",     console_test_cube_quads,        "test.shapes.cube_q"     );
    me->register_console_function( "test.shapes.basic_PAM",  console_test_basic_PAM,         "test.shapes.basic_PAM"  );
    me->register_console_function( "test.build_patches",     console_build_patches,          "test.build_patches"     );
    me->register_console_function( "test.save_patched_obj",  console_save_patched_obj,       "test.save_patched_obj"  );
}
        
void register_test_console_funcs( GLGraphics::MeshEditor* me )
{
    me->register_console_function( "test.set_host_and_module",   console_test_set_host_and_module, "test.set_host_and_module"     );
    me->register_console_function( "test.randomly_rotate",      console_test_randomly_rotate, "test.randomly_rotate"        );
    me->register_console_function( "test.translate_outside_bsphere",   console_test_translate_outside_bsphere, "test.translate_outside_bsphere"        );
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
    
    me->register_console_function(
      "test.save_to_results_folder",
      console_test_save_to_results_folder,
      "test.save_to_results_folder"        );
}

}}