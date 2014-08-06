//
//  PM_basic_console_functions.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 06/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "PM_basic_console_functions.h"
#include <GEL/GLGraphics/MeshEditor.h>
//#include <strstream>
//#include <istream>
//#include <fstream>
//#include <string>
//#include <regex>
//#include <set>
//#include <queue>
//#include "Test.h"
//#include "polarize.h"
#include <MeshEditE/Procedural/Operations/basic_shapes.h>

using namespace GLGraphics;
using namespace Procedural::Operations;

using namespace std;
using namespace HMesh;

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


namespace Procedural
{
namespace ConsoleFuncs
{

    void register_basic_console_funcs(GLGraphics::MeshEditor* me)
    {
        me->register_console_function( "shapes.cube_t",     console_test_cube_triangles,    "shapes.cube_t"     );
        me->register_console_function( "shapes.cube_q",     console_test_cube_quads,        "shapes.cube_q"     );
        me->register_console_function( "shapes.basic_PAM",  console_test_basic_PAM,         "shapes.basic_PAM"  );

    }

}
}