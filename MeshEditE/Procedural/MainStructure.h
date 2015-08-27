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

#include <GEL/HMesh/Manifold.h>
#include <stdio.h>
#include <set>

//#include <MeshEditE/Procedural/Matches/graph_match.h>
#include <GEL/HMesh/Manifold.h>

namespace Procedural {
    
typedef int ModuleID;
typedef std::map< HMesh::VertexID, ModuleID>    VertexModuleMap;

struct StructurePoleInfo{
    PoleInfo    pi;            // pos and normal should be modified here
    ModuleID    moduleID;      // id of the module inside the structure
};
    
struct GluedModuleInfo{
    Module          *module;
    size_t          t_start;
    size_t          connection_valence;
};

class MainStructure{
    
public:
    MainStructure();
    // those methods work only on a logical basis not on a geometrical one
    void glueModule( const HMesh::Manifold& mani,  Module &m, std::vector<Match> &matches  );
    void reAlignIDs( HMesh::VertexIDRemap &remapper );
    void setBoundingSphere( CGLA::Vec3d center, double radius );
    bool isColliding( const Module& m ) const;
    void saveSkeleton( std::string path ) const;
    void saveBVH( std::string path ) const;
    void updatePoleInfo( const HMesh::Manifold& mani );
    
    const Procedural::PoleList& getPoles() const;
    const Procedural::PoleList& getFreePoles() const;
    const Procedural::PoleList& getGluedPoles() const;
    const Procedural::PoleSet&  getFreePoleSet() const;
    const PoleInfo&             getPoleInfo( HMesh::VertexID p ) const;
    inline const PoleInfoMap&   getPoleInfoMap() const{ return freePoleInfoMap;}
    
private:
/************************************************
 * ATTRIBUTES                                   *
 ***********************************************/
    std::vector< GluedModuleInfo >  modules;
    Procedural::PoleList            freePoles;
    Procedural::PoleList            gluedPoles;
    PoleSet                         freePolesSet;
    size_t                          time;
    PoleInfoMap                     freePoleInfoMap;
    Skeleton                        *skel;

};

}

#endif /* defined(__MeshEditE__MainStructure__) */
