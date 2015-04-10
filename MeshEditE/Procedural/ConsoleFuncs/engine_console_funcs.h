//
//  engine_console_funcs.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__engine_console_funcs__
#define __MeshEditE__engine_console_funcs__

#include <stdio.h>

namespace GLGraphics {
    class MeshEditor;
}

namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_engine_console_funcs(GLGraphics::MeshEditor* me);
}}


#endif /* defined(__MeshEditE__engine_console_funcs__) */
