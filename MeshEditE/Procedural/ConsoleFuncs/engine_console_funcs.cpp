//
//  engine_console_funcs.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "engine_console_funcs.h"

#include <sstream>

#include <GEL/GLGraphics/MeshEditor.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/obj_load.h>

#include <MeshEditE/Procedural/StatefulEngine.h>
#include <MeshEditE/Procedural/Toolbox.h>



using namespace std;

using namespace HMesh;
using namespace GLGraphics;

using namespace Procedural::Engines;

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
        a0 >> no_glueings;
    }
    
    if( args.size() > 1 ){
        istringstream a1( args[1] );
        a1 >> no_tests;
    }

    if( no_tests < 1 ) { no_tests = BASE_NO_TESTS; }
    if( no_glueings < 1 ) { no_glueings = BASE_NO_GLUEINGS; }
    
    cout << "testing " << no_tests << " times " << " for " << no_glueings << " glueings" << endl;
    
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.testMultipleTransformations( no_tests, no_glueings );
//    s.alignUsingBestMatch();
}

void align_module( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    //    s.glueModuleToHost();
    s.alignUsingBestMatch();
    
}


void glue_current( MeshEditor *me, const std::vector< std::string > &args ){
        StatefulEngine &s = StatefulEngine::getCurrentEngine();
//    s.glueModuleToHost();
    s.actualGlueing();
    
}


void load_toolbox( MeshEditor *me, const std::vector< std::string > &args ){
    stringstream oss;
    string folder   = "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/Toolbox/";
    string filename = "toolbox1.json";
    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> filename;
    }
    
    string full_path = folder + filename;
    std::cout << "path : "  << full_path;
    
    Procedural::Toolbox& t = Procedural::Toolbox::getToolboxInstance();
    t.clear();
    t.fromJson( full_path );
}



/************************************************
 * DEBUG CALLS                                  *
 ***********************************************/

void art( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.applyRandomTransform();

}

void apa( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.applyOptimalAlignment();

}

void amnth( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.alignModuleNormalsToHost();

}

namespace Procedural{
    namespace ConsoleFuncs{

void register_engine_console_funcs( GLGraphics::MeshEditor* me )
{
    me->register_console_function( "engine.set_host",   set_host, "engine.set_host_and_module" );
    me->register_console_function( "engine.set_module",   set_module, "engine.set_host_and_module" );
    me->register_console_function( "engine.optimal_transform", optimal_transform, "engine.optimal_transform" );
    me->register_console_function( "engine.align", align_module, "engine.align" );
    me->register_console_function( "engine.glue_current", glue_current, "engine.glue_current" );
    me->register_console_function( "engine.load_toolbox", load_toolbox, "engine.load_toolbox" );

    
    // experimental - debug purposes
    me->register_console_function( "engine.debug.apply_random_transform", art, "engine.debug.apply_random_transform" );
    me->register_console_function( "engine.debug.apply_optimal_alignment", apa, "engine.debug.apply_optimal_alignment" );
    me->register_console_function( "engine.debug.align_Module_normals_to_host", amnth, "engine.debug.align_Module_normals_to_host" );
    
}
}}