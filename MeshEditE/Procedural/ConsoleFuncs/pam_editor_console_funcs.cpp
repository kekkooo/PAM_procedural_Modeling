//
//  pam_editor_console_funcs.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 12/02/16.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#include "pam_editor_console_funcs.hpp"

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
#include <MeshEditE/Procedural/pam_editor.hpp>

#include<MeshEditE/Procedural/Helpers/geometric_properties.h>
#include<MeshEditE/Procedural/Helpers/misc.h>
#include<polarize.h>



using namespace std;

using namespace HMesh;
using namespace GLGraphics;
using namespace Util;

using namespace Procedural::Engines;
using namespace Procedural::Geometry;
using namespace Procedural;

#define BASE_NO_TESTS 10
#define BASE_NO_GLUEINGS 1

#define MODULES_FOLDER "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Shapes/Modules/"
#define RESULTS_FOLDER "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Results/"
#define TESTS_FOLDER   "/Users/francescousai/Documents/Dottorato/Conferenze/CGI_PG2015/Tests/"

void load_config( MeshEditor *me, const std::vector< std::string > &args ){
    stringstream oss;
    string folder   = MODULES_FOLDER;
    string filename = "toolbox1.json";
    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> filename;
    }
    
    string full_path = folder + filename;
    std::cout << "path : "  << full_path;
    
    Procedural::PamEditor& p = Procedural::PamEditor::getCurrentPamEditor();
    p.set_config( full_path );
    
}


namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_pam_editor_console_funcs( GLGraphics::MeshEditor* me )
        {
            me->register_console_function( "pam_editor.load_config", load_config, "pam_editor.load_config" );
        }
    }}