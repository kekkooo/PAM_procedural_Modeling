//
//  engine_console_funcs.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "engine_console_funcs.h"

#include <sstream>
//#include <random>

#include <GEL/GLGraphics/MeshEditor.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/obj_load.h>

#include <MeshEditE/Procedural/StatefulEngine.h>

//#include <MeshEditE/Procedural/Operations/basic_shapes.h>
//#include <MeshEditE/Procedural/Operations/Algorithms.h>
//#include <MeshEditE/Procedural/Operations/structural_operations.h>
//#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
//#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
//#include <MeshEditE/Procedural/PMEngine.h>
//#include <MeshEditE/Procedural/Matches/Matches.h>
//#include <MeshEditE/Procedural/Helpers/manifold_copy.h>

//#include <MeshEditE/Test.h>
//#include <MeshEditE/Tests/test_random_rotation.h>

//#include "patch_mapping.h"
//#include "polarize.h"

using namespace std;

using namespace HMesh;
using namespace GLGraphics;

using namespace Procedural::Engines;
//using namespace Procedural::Operations;
//using namespace Procedural::Operations::Algorithms;

#define BASE_NO_TESTS 10
#define BASE_NO_GLUEINGS 1

#define MODULES_FOLDER "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/Modules/"

void set_host( MeshEditor *me, const std::vector< std::string > &args ){
    
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.setHost( me->active_mesh( ));
}


void set_module( MeshEditor *me, const std::vector< std::string > &args ){

    Manifold module;
    stringstream oss;
    string filename = "1.obj";
    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> filename;
    }
    
    oss << MODULES_FOLDER << filename;
    obj_load( oss.str(), module );
    
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.setModule( module );
}

void optimal_transform( MeshEditor *me, const std::vector< std::string > &args ){
    
    stringstream oss;
    int no_tests = BASE_NO_TESTS;
    int no_glueings = BASE_NO_GLUEINGS;

    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> no_tests;
    }
    
    if( args.size() > 1 ){
        istringstream a1( args[1] );
        a1 >> no_glueings;
    }

    if( no_tests < 1 ) { no_tests = BASE_NO_TESTS; }
    if( no_glueings < 1 ) { no_glueings = BASE_NO_GLUEINGS; }
    
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.testMultipleTransformations( no_tests, no_glueings );
    s.alignUsingBestMatch();
}

void glue_current( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
//    s.glueModuleToHost();
    s.actualGlueing();
    
}

void set_constraint_dimensionality( MeshEditor *me, const std::vector< std::string > &args ){
    
    stringstream oss;
    int dim_constraint = 3;
    
    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> dim_constraint;
    }
    
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    
    switch (dim_constraint) {
        case 1:
            s.setConstraint1D();
            break;
        case 2:
            s.setConstraint2D();
            break;
        case 3:
            s.setConstraint3D();
            break;
            
        default:
            me->printf("option not available. only 1( 1D ), 2( 2D ), 3 ( 3D )");
            break;
    }
}

namespace Procedural{
    namespace ConsoleFuncs{

void register_engine_console_funcs( GLGraphics::MeshEditor* me )
{
    me->register_console_function( "engine.set_host",   set_host, "engine.set_host_and_module" );
    me->register_console_function( "engine.set_module",   set_module, "engine.set_host_and_module" );
    me->register_console_function( "engine.optimal_transform", optimal_transform, "engine.optimal_transform" );
    me->register_console_function( "engine.set_constraint_dimensionality", set_constraint_dimensionality, "engine.set_constraint_dimensionality" );
    me->register_console_function( "engine.glue_current", glue_current, "engine.glue_current" );
    
}
}}