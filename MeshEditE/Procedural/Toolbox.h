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
#include <random>

namespace Procedural{
    
    struct ModuleInfo{
        Module*         m           = NULL;
        float           probability = 1.0;
        size_t          no_pieces   = 0;
        size_t          no_glueings = 1;
        std::string     name;
    };
    
class Toolbox{
    
    public :
        static  Toolbox& getToolboxInstance();
    
        bool hasNext()      const;
        Module& getNext();

        void addModule( std::string path, size_t no_pieces );
        void fromJson( std::string path );
        void clear();
        void undoLast();
        void print() const;
    
        inline size_t noRemainingPieces() const { return total_pieces; }

private :
    Toolbox();
    Toolbox( Toolbox const& ) = delete;
    void operator   = (Toolbox const&)               = delete;
    void updateProbabilities();
    
    
private :

    std::vector<ModuleInfo*> modules;
    size_t                  total_pieces;
    std::mt19937_64         randomizer;
    float                   rand_max;
    size_t                  last_used_module;
    bool                    used_module = false;

//      std::function<int(float)> random_to_module;

};
}
#endif /* defined(__MeshEditE__Toolbox__) */
