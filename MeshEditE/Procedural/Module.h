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
#include <queue>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>

#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Mat4x4d.h>

#include <polarize.h>

#include "pam_skeleton.h"
#include "collision_detection.h"

namespace Procedural{
    
typedef std::pair< HMesh::VertexID, HMesh::VertexID >   Match;
    
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
    
struct PoleAnisotropyInfo{
    HMesh::VertexID directionID  = HMesh::InvalidVertexID;
    CGLA::Vec3d     direction;
    bool            is_defined   = false;
    bool            is_bilateral = false ;
    
    std::string toString() const{
        std::stringstream oss;
        if( is_defined ){
        oss << "\t ID  : " << directionID << std::endl
            << "\t dir : " << direction << std::endl
            << "\t Bi  : " << (is_bilateral ? " yes " : "no") << std::endl;
        }
        else{ oss << "not defined" << std::endl; }
        return oss.str();
    }
};
    
struct PoleGeometryInfo{
    unsigned int    valence;
    CGLA::Vec3d     pos;
    CGLA::Vec3d     normal;

    std::string toString() const{
        std::stringstream oss;
        oss << "\t val   : " << valence << std::endl
            << "\t pos   : " << pos << std::endl
            << "\t norm  : " << normal << std::endl;
        return oss.str();
    }
};
    

    
struct PoleInfo{
    HMesh::VertexID     original_id;
    PoleAnisotropyInfo  anisotropy;
    PoleGeometryInfo    geometry;
    Moduletype          moduleType  = 0;     // an identifier for the module's type

    int                 age;                 // pole's starting age
    bool                isFree      = true;
    bool                isActive    = true;
    bool                can_connect_to_self = true; /* can connect to an instance of the 
                                                     same pole on another module of the 
                                                     same exact type ( same file ) */
    
    std::string toString() const{
        std::stringstream oss;
        
        oss
        << " moduletype : " << moduleType << std::endl
        << "original ID : " << original_id << std::endl
        << "anisotropy  : "
        << " [ " << std::endl << anisotropy.toString() << " ] " << std::endl
        << "geometry [ " << std::endl << geometry.toString() << " ] " << std::endl
        << "active      : " << (isActive ? " yes " : "no") << std::endl
        << "self-connect: " << (can_connect_to_self ? " yes " : "no") << std::endl;

        return oss.str();
    }
};
    
typedef std::map<HMesh::VertexID, PoleInfo>             PoleInfoMap;
typedef std::vector< HMesh::VertexID>                   PoleList;
typedef std::set< HMesh::VertexID >                     PoleSet;

class Module{
    
public:
         // this will instantiate the internal manifold structure and pole info using obj_load
        Module () { m = NULL; }
        Module( std::string path, std::string config, Moduletype mType );
        Module( HMesh::Manifold &manifold, Moduletype mType );
    
        Module& getTransformedModule( const CGLA::Mat4x4d &T, bool transform_geometry = false );
        void reAlignIDs( HMesh::VertexIDRemap &remapper );
    
        const PoleInfo&    getPoleInfo( HMesh::VertexID p ) const;
        bool isPole( HMesh::VertexID v );
        const Skeleton& getSkeleton() const;

        inline const PoleInfoMap& getPoleInfoMap()const{ return poleInfoMap; }
    
        static bool poleCanMatch( const PoleInfo& p1, const PoleInfo& p2);
    
    
        void sanityCheck();
    
private:
    void    BuildPoleInfo();
    void    LoadPoleConfig( std::string path );
    void    getPoleAnisotropy( HMesh::VertexID pole, CGLA::Vec3d& dir,  HMesh::VertexID neighbor ) const;
    
    
/************************************************
* ATTRIBUTES                                   *
***********************************************/
public:
#warning those should be deleted from public
    HMesh::Manifold     *m;
    int                 no_of_glueings;
    PoleList            poleList;
    PoleSet             poleSet;

    CGLA::Vec3d         bsphere_center;
    double              bsphere_radius;
    
private :
    PoleInfoMap         poleInfoMap;
    Skeleton            *skeleton;

    };
}

#endif /* defined(__MeshEditE__Module__) */
