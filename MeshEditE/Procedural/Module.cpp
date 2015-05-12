//
//  Module.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "Module.h"
#include "Test.h"

using namespace HMesh;

namespace Procedural{

Module::Module( std::string path, Moduletype mType ){
    this->m = NULL;
}
    
Module Module::getTransformedModule( const CGLA::Mat4x4d &T )
{
    Module M;// = new Module();
    
    M.no_of_glueings = this->no_of_glueings;
    M.bsphere_center = T.mul_3D_point( bsphere_center );
    M.bsphere_radius = bsphere_radius;
    
    for( VertexID vid : this->poleList ){
        M.poleList.push_back( vid );
        M.poleInfoMap[vid].geometry.pos     = T.mul_3D_point( poleInfoMap[vid].geometry.pos );
        M.poleInfoMap[vid].geometry.normal  = mul_3D_dir( T, poleInfoMap[vid].geometry.normal );
        M.poleInfoMap[vid].geometry.valence = poleInfoMap[vid].geometry.valence;
        M.poleInfoMap[vid].labels           = poleInfoMap[vid].labels;
        M.poleInfoMap[vid].age              = poleInfoMap[vid].age;
        M.poleInfoMap[vid].isFree           = poleInfoMap[vid].isFree;
    }
    
    assert( M.poleInfoMap.size() == this->poleInfoMap.size( ));
    
    return M;
}
    

}