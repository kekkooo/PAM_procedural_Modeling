//
//  Module.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "Module.h"

using namespace HMesh;

namespace Procedural{

Module::Module( std::string path, Moduletype mType ){
    this->m = NULL;
}
    
void Module::applyTransformToPoleInfo( CGLA::Mat4x4d &T )
{
    for( auto it = poleInfoMap.begin(); it != poleInfoMap.end(); ++it )
    {
        VertexID key = it->first;
        poleInfoMap[key].geometry.pos       = T.mul_3D_point( poleInfoMap[key].geometry.pos );
        poleInfoMap[key].geometry.normal    = T.mul_3D_point( poleInfoMap[key].geometry.normal );
    }
}
    

}