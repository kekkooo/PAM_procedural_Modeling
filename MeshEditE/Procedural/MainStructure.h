//
//  MainStructure.h
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__MainStructure__
#define __MeshEditE__MainStructure__

#include "Module.h"

#include <stdio.h>
#include <set>

#include <MeshEditE/Procedural/Matches/graph_match.h>
#include <GEL/HMesh/Manifold.h>

namespace Procedural {
    
typedef int ModuleID;
typedef std::map< HMesh::VertexID, ModuleID>    VertexModuleMap;
typedef std::set< HMesh::VertexID >             PoleSet;

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
    void reAlignIDs( HMesh::VertexIDRemap &remapper );
    const Procedural::PoleList& getPoles();
    const Procedural::PoleList& getFreePoles();
    const Procedural::PoleList& getGluedPoles();
    const PoleSet&              getFreePoleSet();

    
    
    
private:
/************************************************
 * ATTRIBUTES                                   *
 ***********************************************/
    std::vector< GluedModuleInfo >  modules;
    Procedural::PoleList            freePoles;
    Procedural::PoleList            gluedPoles;
    PoleSet                         freePolesSet;
    size_t                          time;

};

}

#endif /* defined(__MeshEditE__MainStructure__) */
