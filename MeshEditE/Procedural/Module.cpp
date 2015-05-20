//
//  Module.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "Module.h"
#include "Test.h"
#include "MeshEditE/Procedural/Helpers/geometric_properties.h"
#include <GEL/HMesh/obj_load.h>
#include <fstream>

using namespace HMesh;
using namespace CGLA;
using namespace Procedural::Geometry;
using namespace std;

namespace Procedural{

Module::Module( std::string path, Moduletype mType ){
    cout << "trying to load file : " << path << endl;
    ifstream f( path );
    assert( f.good() );
    
    this->m = new Manifold();
    
    obj_load( path, *this->m );
    bsphere( *m, bsphere_center, bsphere_radius );
    BuildPoleInfo();
}
    
void Module::BuildPoleInfo(){
    
    assert( m != NULL );
    
    for( VertexID vid : m->vertices() ){
        if( is_pole( *m, vid )){
            PoleInfo pi;
            pi.geometry.valence = valence( *m, vid );
            pi.geometry.pos     = m->pos( vid );
            Vec3d n = vertex_normal( *m, vid );
            n.normalize();
            pi.geometry.normal = n;
            this->poleList.push_back( vid );
            this->poleInfoMap[vid] = pi;
        }
    }
}
    
    
Module::Module( Manifold &manifold, Moduletype mType, size_t no_glueings ){
    this->m = &manifold;
    Vec3d centroid;
    double radius;
    bsphere( *this->m, centroid, radius );
    this->bsphere_center = centroid;
    this->bsphere_radius = radius;
    this->no_of_glueings = no_glueings;
    
    BuildPoleInfo();
}
    
Module& Module::getTransformedModule( const CGLA::Mat4x4d &T, bool transform_module )
{
    Module *M = new Module();
    M->m             = this->m;
    M->poleList      = PoleList();
    M->poleInfoMap   = PoleInfoMap();
    
    M->no_of_glueings = this->no_of_glueings;
    M->bsphere_center = T.mul_3D_point( bsphere_center );
    M->bsphere_radius = bsphere_radius;
    
    for( VertexID vid : this->poleList ){
        M->poleList.push_back( vid );
        M->poleInfoMap[vid].geometry.pos     = T.mul_3D_point( poleInfoMap[vid].geometry.pos );
        M->poleInfoMap[vid].geometry.normal  = mul_3D_dir( T, poleInfoMap[vid].geometry.normal );
        M->poleInfoMap[vid].geometry.valence = poleInfoMap[vid].geometry.valence;
        M->poleInfoMap[vid].labels           = poleInfoMap[vid].labels;
        M->poleInfoMap[vid].age              = poleInfoMap[vid].age;
        M->poleInfoMap[vid].isFree           = poleInfoMap[vid].isFree;
    }
    
    if( transform_module ){
        for( VertexID v : m->vertices()){
            m->pos( v ) = T.mul_3D_point( m->pos( v ));
        }
    }
    
    assert( M->poleInfoMap.size() == this->poleInfoMap.size( ));
    
    return *M;
}
    
    void Module::reAlignIDs(HMesh::VertexIDRemap &remapper){
        // realign poleList
        PoleInfoMap p;
        for( int i = 0; i < poleList.size(); ++i ){
            VertexID newID = remapper[poleList[i]];
            p[newID] = poleInfoMap[poleList[i]];
            poleList[i] = newID;
        }
        poleInfoMap = std::move( p );
        
    }
    

}