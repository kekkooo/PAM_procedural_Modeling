//
//  engine_console_funcs.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "engine_console_funcs.h"

#include <unistd.h>
#include <sstream>
#include <chrono>

#include <GEL/GLGraphics/MeshEditor.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/obj_load.h>
#include <GEL/Util/Timer.h>

#include <MeshEditE/Procedural/StatefulEngine.h>
#include <MeshEditE/Procedural/Toolbox.h>
#include <MeshEditE/Procedural/Helpers/manifold_copy.h>

#include<MeshEditE/Procedural/Helpers/misc.h>


using namespace std;

using namespace HMesh;
using namespace GLGraphics;
using namespace Util;

using namespace Procedural::Engines;
using namespace Procedural;

#define BASE_NO_TESTS 10
#define BASE_NO_GLUEINGS 1

#define MODULES_FOLDER "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/Modules/"
#define RESULTS_FOLDER "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Results/"
#define TESTS_FOLDER   "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Tests/"

static size_t               current_conf = 0;
std::vector<CGLA::Mat4x4d>  ts;
CGLA::Mat4x4d               current_T;
VertexSet                   M_vertices;
std::string                 curr_toolbox;


long millis(){
    return chrono::system_clock::now().time_since_epoch().count();;
}

void set_host( MeshEditor *me, const std::vector< std::string > &args ){
    
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.setHost( me->active_mesh( ));
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
    s.testMultipleTransformations( );
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
        curr_toolbox = filename;
    }
    
    string full_path = folder + filename;
    std::cout << "path : "  << full_path;
    
    Procedural::Toolbox& t = Procedural::Toolbox::getToolboxInstance();
    t.clear();
    t.fromJson( full_path );
}

struct toolbox_step_result{
    bool has_next           = true;
    bool can_glue           = true;
    bool enough_free_poles  = true;
    inline bool ok() { return has_next && can_glue && enough_free_poles ; }
};

toolbox_step_result toolbox_step( Procedural::Toolbox &t, StatefulEngine &s ){
    toolbox_step_result result;
    result.has_next = t.hasNext();
    if( !result.ok() ){ return result; }

    Timer timer;
    cout << endl << "########################################" << endl << endl;
    timer.start();
    Module m = t.getNext();
    float t0 = timer.get_secs();
    cout << " adding a piece with "  << m.no_of_glueings << "-valent connection in " << t0 << "s" << endl;
    
    s.setModule( m );
    
    float t1 = timer.get_secs();
    
    cout << "module set in : " << (t1-t0) << endl;
    
    result.enough_free_poles   = s.noFreePoles() >= m.no_of_glueings;
    if( result.enough_free_poles){
        float t1_1 = timer.get_secs();
        result.can_glue  = s.testMultipleTransformations();
        float t2 = timer.get_secs();
        cout << "configurations tested in in : " << (t2-t1_1) << "s" << endl;
        
    }
    
    if( result.can_glue && result.enough_free_poles ){
        float t2_1 = timer.get_secs();
        s.glueCurrent();
        float t3 = timer.get_secs();
        cout << "glueing done in : " << (t3-t2_1) << "s" << endl;
    }
    else{
        t.undoLast();
    }
    return result;
}

void toolbox_handle_result( const Procedural::Toolbox &t, const toolbox_step_result& result ){
    if( !result.has_next ){
        cout << "no more pieces " << endl;
    }
    if( !result.can_glue ){
        cout << "cannot find a feasible solution. remaining pieces : " << t.noRemainingPieces() << endl;
    }
    if( !result.enough_free_poles ){
        cout << "there aren't enough free poles. remaining pieces :" << t.noRemainingPieces() << endl;
    }
    t.print();

}

void empty_toolbox( MeshEditor *me, const std::vector< std::string > &args ){
    
    Procedural::Toolbox& t = Procedural::Toolbox::getToolboxInstance();
    StatefulEngine &s = StatefulEngine::getCurrentEngine();

    toolbox_step_result result;
    
    while( result.ok() ){
        result = toolbox_step( t, s );
    }
    toolbox_handle_result( t, result );
}


void step_toolbox( MeshEditor *me, const std::vector< std::string > &args ){
    
    Procedural::Toolbox& t = Procedural::Toolbox::getToolboxInstance();
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    
    toolbox_step_result result = toolbox_step( t, s );
    if( !result.has_next ){
        cout << "no more pieces " << endl;
    }
    if( !result.can_glue ){
        cout << "cannot find a feasible solution. remaining pieces : " << t.noRemainingPieces() << endl;
    }
    if( !result.enough_free_poles ){
        cout << "there aren't enough free poles. remaining pieces :" << t.noRemainingPieces() << endl;
    }
    toolbox_handle_result( t, result );
}




/************************************************
 * DEBUG CALLS                                  *
 ***********************************************/

void clean_conf( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    Procedural::Helpers::test_delete( *(s.m), M_vertices );
}

void gen_confs( MeshEditor *me, const std::vector< std::string > &args ){
    Procedural::Toolbox& t = Procedural::Toolbox::getToolboxInstance();
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    
    s.setModule( t.getNext() );

    ts.clear();
    s.buildTransformationList( ts );
    current_conf = 0;
    assert(s.transformedModules.size() > 0 );
}

void save_all_configurations( MeshEditor *me, const std::vector< std::string > &args ){
    Procedural::Toolbox& t = Procedural::Toolbox::getToolboxInstance();
    StatefulEngine &s = StatefulEngine::getCurrentEngine();    
    s.setModule( t.getNext() );
    Procedural::Module *m = s.candidateModule;
    
    long ms = millis();

    ts.clear();
    s.buildTransformationList( ts );
    
    if( ts.size() <= 0 ){ cout << "there are no transformations available"; return; }
    stringstream oss;
    oss << curr_toolbox << '_' << ms;
    string folder_name = oss.str();
    string full_path = TESTS_FOLDER + oss.str();
    Procedural::Helpers::Misc::new_folder( TESTS_FOLDER, folder_name );
    
    
    size_t count = 0;
    for( const auto& T : ts ){
        Manifold *tm = new Manifold();
        VertexIDRemap ___, ____;
        
        // copy geometry to clean manifold
        Procedural::Helpers::add_manifold( *tm, *(m->m), ___, ___, M_vertices );
        // apply T
        for( VertexID vid :  tm->vertices() ){
            tm->pos( vid ) = T.mul_3D_point(tm->pos(vid));
        }
        
        ___.clear();
        ____.clear();
        M_vertices.clear();
        // copy transformed geometry to host
        Procedural::Helpers::add_manifold( *(s.m), *(tm), ___, ___, M_vertices );
        //save
        stringstream oss;
        oss << full_path << "/" << count << ".obj";
        obj_save( oss.str(), *s.m );
//        usleep(1000000);

        //remove copied geometry
        Procedural::Helpers::test_delete( *(s.m), M_vertices );
        ++count;
    }
    
    
}

void next_configuration( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    Procedural::Module *m = s.candidateModule;
    
    if( M_vertices.size() > 0 ){
            Procedural::Helpers::test_delete( *(s.m), M_vertices );
    }
    
    
    if( current_conf >= ts.size()) { cout << "end reached" << endl; return; }

    CGLA::Mat4x4d t = ts[++current_conf];
    
    Manifold *tm = new Manifold();

    VertexIDRemap ___, ____;
    
    Procedural::Helpers::add_manifold( *tm, *(m->m), ___, ___, M_vertices );
    for( VertexID vid :  tm->vertices() ){
        tm->pos( vid ) = t.mul_3D_point(tm->pos(vid));
    }
    
    ___.clear();
    ____.clear();
    M_vertices.clear();
    
    Procedural::Helpers::add_manifold( *(s.m), *(tm), ___, ___, M_vertices );

}

void prev_configuration( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    Procedural::Module *m = s.candidateModule;

    CGLA::Mat4x4d t = ts[--current_conf];
    Procedural::Module tm = m->getTransformedModule( t, true );
    
    VertexIDRemap ___, ____;
    M_vertices.clear();
    
    Procedural::Helpers::add_manifold( *(s.m), *(tm.m), ___, ___, M_vertices );
}

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


void saveSkel( MeshEditor *me, const std::vector< std::string > &args ){
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.getMainStructure().saveSkeleton("//Users//francescousai//Desktop//example.skel");
}

void saveBVH( MeshEditor *me, const std::vector< std::string > &args ){
    //    "//Users//francescousai//Desktop//bvh.skel"
    StatefulEngine &s = StatefulEngine::getCurrentEngine();
    s.getMainStructure().saveBVH("//Users//francescousai//Desktop//bvh.skel");
}

namespace Procedural{
    namespace ConsoleFuncs{
        

void register_engine_console_funcs( GLGraphics::MeshEditor* me )
{
    me->register_console_function( "engine.set_host",   set_host, "engine.set_host_and_module" );
    me->register_console_function( "engine.optimal_transform", optimal_transform, "engine.optimal_transform" );
    me->register_console_function( "engine.align", align_module, "engine.align" );
    me->register_console_function( "engine.glue_current", glue_current, "engine.glue_current" );

    me->register_console_function( "engine.toolbox.load", load_toolbox, "engine.toolbox.load" );
    me->register_console_function( "engine.toolbox.empty", empty_toolbox, "engine.toolbox.empty" );
    me->register_console_function( "engine.toolbox.step", step_toolbox, "engine.toolbox.step" );

    
    // experimental - debug purposes
    
    me->register_console_function( "engine.debug.next_configuration", next_configuration, "engine.debug.next_configuration" );
    me->register_console_function( "engine.debug.prev_configuration", prev_configuration, "engine.debug.prev_configuration" );
    me->register_console_function( "engine.debug.gen_conf", gen_confs, "engine.debug.gen_conf" );
    me->register_console_function( "engine.debug.clean_conf", clean_conf, "engine.debug.clean_conf" );
    me->register_console_function( "engine.debug.save_confs", save_all_configurations, "engine.debug.save_confs" );
    
    me->register_console_function( "engine.debug.apply_random_transform", art, "engine.debug.apply_random_transform" );
    me->register_console_function( "engine.debug.apply_optimal_alignment", apa, "engine.debug.apply_optimal_alignment" );
    me->register_console_function( "engine.debug.align_Module_normals_to_host", amnth, "engine.debug.align_Module_normals_to_host" );
    me->register_console_function( "engine.debug.saveSkel", saveSkel, "engine.debug.saveSkel" );
    me->register_console_function( "engine.debug.saveBVH", saveBVH, "engine.debug.saveBVH" );
    
}
}}