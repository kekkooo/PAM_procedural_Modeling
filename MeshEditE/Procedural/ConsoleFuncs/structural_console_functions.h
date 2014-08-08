//
//  structural_console_functions.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__structural_console_functions__
#define __MeshEditE__structural_console_functions__

#include <iostream>

namespace GLGraphics {
    class MeshEditor;
}

namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_structural_console_funcs(GLGraphics::MeshEditor* me);
        
}}


#endif /* defined(__MeshEditE__structural_console_functions__) */
