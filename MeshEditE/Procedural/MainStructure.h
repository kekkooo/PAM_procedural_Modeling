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

namespace Procedural {
    
typedef int ModuleID;
typedef std::map< HMesh::VertexID, ModuleID> VertexModuleMap;

struct StructurePoleInfo{
    PoleInfo    pi;            // pos and normal should be modified here
    ModuleID    moduleID;      // id of the module inside the structure
};
    
struct GluedModuleInfo{
    Module          *module;
    int             t_start;
    CGLA::Vec3d     centroid;
    double          radius;
    int             connection_valence;
};

class MainStructure{
    
public:
    MainStructure() {};
    
private:
/************************************************
 * ATTRIBUTES                                   *
 ***********************************************/
    std::vector< GluedModuleInfo >  modules;

};

}

#endif /* defined(__MeshEditE__MainStructure__) */
