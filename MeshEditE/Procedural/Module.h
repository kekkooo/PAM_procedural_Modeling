//
//  Module.h
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Module__
#define __MeshEditE__Module__

#include <stdio.h>

#include <map>
#include <vector>

#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>

namespace Procedural{
    
enum PoleLabel {
    NoneP    = 0,       // this is 00000000 it will match no value
    RedP     = 1,
    BlueP    = 2,
    GreenP   = 4,
    CyanP    = 8,
    MagentaP = 16,
    YellowP  = 32,
    PurpleP  = 64,
    AllP     = 127      // this is 11111111 it will match each value
};

typedef unsigned int Moduletype;
    
struct PoleGeometryInfo{
    unsigned int    valence;
    CGLA::Vec3d     pos;
    CGLA::Vec3d     normal;
};
    
struct PoleInfo{
    PoleGeometryInfo    geometry;
    Moduletype          moduleType  = 0;     // an identifier for the module's type
    PoleLabel           labels      = AllP;  // a label atteched to a pole
    int                 age;                 // pole's starting age
    bool                isFree      = false;
};
    
typedef std::map<HMesh::VertexID, PoleInfo>             PoleInfoMap;



class Module{
    
public:
    Module( std::string path, Moduletype mType ); // this will instantiate the internal manifold structure and pole info using obj_load
#warning those should be deleted from public
    PoleInfoMap         poleInfoMap;
    CGLA::Vec3d         bsphere_center;
    double              bsphere_radius;

    
private:
    void BuildPoleInfo();
    
        
/************************************************
* ATTRIBUTES                                   *
***********************************************/
private:
    HMesh::Manifold     *m;
    //PoleInfoMap         poleInfoMap;
    int                 no_of_glueings;
//    CGLA::Vec3d         bsphere_center;
//    double              bsphere_radius;

    
        
    };
}

#endif /* defined(__MeshEditE__Module__) */