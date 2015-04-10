//
//  PM_basic_console_functions.h
//  MeshEditE
//
//  Created by Francesco Usai on 06/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__PM_basic_console_functions__
#define __MeshEditE__PM_basic_console_functions__

#include <iostream>

namespace GLGraphics {
    class MeshEditor;
}

namespace Procedural{
    namespace ConsoleFuncs{
    
    void register_basic_console_funcs(GLGraphics::MeshEditor* me);
    void register_algorithm_console_funcs(GLGraphics::MeshEditor* me);
    void register_match_console_funcs(GLGraphics::MeshEditor* me);
        void register_test_console_funcs( GLGraphics::MeshEditor* me );
}}

#endif /* defined(__MeshEditE__PM_basic_console_functions__) */
