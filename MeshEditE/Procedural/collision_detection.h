//
//  collision_detection.h
//  MeshEditE
//
//  Created by Francesco Usai on 12/06/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef MeshEditE_collision_detection_h
#define MeshEditE_collision_detection_h

#include "pam_skeleton.h"

#include <set>

namespace Procedural{
    
    inline bool ball_collide( const Ball& b1, const Ball& b2 ){
        double dist = ( b1.center - b2.center ).length();
        return ( dist < b1.radius + b2.radius );
    }
    
    bool collide( const Skeleton& main, const Skeleton& other );
}

#endif