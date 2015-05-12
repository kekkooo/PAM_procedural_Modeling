//
//  MainStructure.h
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__MainStructure__
#define __MeshEditE__MainStructure__

#include <stdio.h>

#include "Module.h"
#include <MeshEditE/Procedural/Matches/graph_match.h>

namespace Procedural {
    
typedef int ModuleID;
typedef std::map< HMesh::VertexID, ModuleID> VertexModuleMap;

struct StructurePoleInfo{
    PoleInfo    pi;            // pos and normal should be modified here
    ModuleID    moduleID;      // id of the module inside the structure
};
    
struct GluedModuleInfo{
    Module          *module;
    size_t          t_start;
    int             connection_valence;
};

class MainStructure{
    
public:
    MainStructure();
    // those methods work only on a logical basis not on a geometrical one
    void glueModule( Module &m, std::vector<Procedural::GraphMatch::Match> &matches  );
    const Procedural::PoleList& getPoleList();
    
    
    
private:
/************************************************
 * ATTRIBUTES                                   *
 ***********************************************/
    std::vector< GluedModuleInfo >  modules;
    Procedural::PoleList            freePoles;
    Procedural::PoleList            gluedPoles;
    size_t                          time;

};

}

#endif /* defined(__MeshEditE__MainStructure__) */
