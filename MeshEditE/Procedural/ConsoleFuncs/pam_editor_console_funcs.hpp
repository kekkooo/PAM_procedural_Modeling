//
//  pam_editor_console_funcs.hpp
//  MeshEditE
//
//  Created by Francesco Usai on 12/02/16.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef pam_editor_console_funcs_hpp
#define pam_editor_console_funcs_hpp

#include <stdio.h>

namespace GLGraphics {
    class MeshEditor;
}

namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_pam_editor_console_funcs(GLGraphics::MeshEditor* me);
    }}


#endif /* pam_editor_console_funcs_hpp */
