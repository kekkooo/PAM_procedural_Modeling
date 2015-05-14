//
//  Toolbox.h
//  MeshEditE
//
//  Created by Francesco Usai on 14/05/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Toolbox__
#define __MeshEditE__Toolbox__

#include <stdio.h>
#include "Module.h"

namespace Procedural{
    
    struct ModuleInfo{
        Module* m           = NULL;
        float   probability = 1.0;
        size_t  no_pieces   = 0;
        size_t  no_glueings = 1;
    };
    
class Toolbox{
    
    public :
        static  Toolbox& getToolboxInstance();
        void addModule( std::string path, size_t no_pieces );
        const Module& getNext();
        void fromJson( std::string path );
        void clear();

private :
    Toolbox();
    Toolbox( Toolbox const& ) = delete;
    void operator   = (Toolbox const&)               = delete;
    void updateProbabilities();
    
private :

    std::vector<ModuleInfo> modules;
    int                     total_pieces;

};
}
#endif /* defined(__MeshEditE__Toolbox__) */
