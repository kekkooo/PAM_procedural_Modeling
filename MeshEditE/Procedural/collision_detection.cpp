//
//  collision_detection.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 13/06/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include <stdio.h>
#include "collision_detection.h"

namespace Procedural{
bool collide( const Skeleton& main, const Skeleton& other ){
    
    const ShapeBall& main_cd    = main.getCdHierarchy(),
    other_cd   = other.getCdHierarchy();
    
#warning it should be an in-order visit, for now it is a level-visit
    
    // if there is no collision between the bounding spheres I can safely return false
    if( !ball_collide( main_cd.ball, other_cd.ball )){ return false; }
    
    std::vector< std::pair< NodeID, NodeID > >      joints_to_test;
    std::set< std::pair< BoneID, BoneID > >         bones_to_test;
    size_t tests = 0, j_tests = 0, b_tests = 0, n_tests = 0;
    
    
    for( const auto& item : main_cd.joints ){
        for( const auto& other_item : other_cd.joints ){
            ++tests;
            ++j_tests;
            if( ball_collide( item.second.ball, other_item.second.ball )){
                joints_to_test.push_back( std::make_pair( item.first, other_item.first ) );
            }
        }
    }
    // for each pair of joints to test, check which pair of bones should be tested
    for( const auto& item : joints_to_test ){
        for( BoneID bm : main_cd.joints.at( item.first ).incidentBones ){
            assert( main_cd.bones.count( bm ) > 0 );
            const Ball b1 = main_cd.bones.at( bm ).ball;
            
            for( BoneID bo : other_cd.joints.at( item.second ).incidentBones ){
                assert( other_cd.bones.count( bo ) > 0 );
                const Ball b2 = other_cd.bones.at( bo ).ball;
                
                ++tests;
                ++b_tests;
                
                if( ball_collide( b1, b2) ){
                    bones_to_test.insert( std::make_pair( bm, bo ) );
                }
            }
        }
    }
    for( const auto& item : bones_to_test ){
        for( NodeID nm : main_cd.bones.at( item.first ).nodes ){
            assert( nm >= 0 && nm < main.nodes.size( ));
            const Ball b1 = main.nodes.at( nm ).ball;
            
            for( NodeID no : other_cd.bones.at( item.second ).nodes ){
                assert( no >= 0 && no < other.nodes.size( ));
                const Ball b2 = other.nodes.at( no ).ball;
                ++tests;
                ++n_tests;
                if( ball_collide( b1, b2) ){
                    std::cout << "T : " << tests << " # J : " << j_tests <<
                    " # B : " << b_tests << " # N " << n_tests << std::endl;
                    std::cout << " no joints : ( " << main_cd.joints.size()  << ", " << other_cd.joints.size()   << " ) " << std::endl;
                    std::cout << " no bones  : ( " << main_cd.bones.size()   << ", " << other_cd.bones.size()    << " ) " << std::endl;
                    std::cout << " no nodes  : ( " << main.nodes.size()      << ", " << other.nodes.size()       << " ) " << std::endl;
                    
                    return true;
                }
            }
        }
    }
    
    std::cout << "T : " << tests << " # J : " << j_tests <<
    " # B : " << b_tests << " # N " << n_tests << std::endl;
    std::cout << " no joints : ( " << main_cd.joints.size()  << ", " << other_cd.joints.size()   << " ) " << std::endl;
    std::cout << " no bones  : ( " << main_cd.bones.size()   << ", " << other_cd.bones.size()    << " ) " << std::endl;
    std::cout << " no nodes  : ( " << main.nodes.size()      << ", " << other.nodes.size()       << " ) " << std::endl;
    
    
    return false;
}
}
