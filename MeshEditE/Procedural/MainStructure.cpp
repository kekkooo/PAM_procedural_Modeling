#define TRACE
//
//  MainStructure.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "MainStructure.h"

#include <set>
#include <algorithm>

#include "Test.h"

using namespace std;
using namespace HMesh;
using namespace Procedural;
using namespace Procedural::GraphMatch;

namespace Procedural{
    
    MainStructure::MainStructure(){
        time = 0;
        skel = new Skeleton();
    }

    const PoleList& MainStructure::getPoles(){
        return freePoles;
    }

    const PoleList& MainStructure::getFreePoles(){
        return freePoles;
    }
    
    const PoleSet& MainStructure::getFreePoleSet(){
        return freePolesSet;
    }

    const PoleList& MainStructure::getGluedPoles(){
        return freePoles;
    }
    
    const PoleInfo& MainStructure::getPoleInfo( HMesh::VertexID p ){
        assert(freePoleInfoMap.count(p) > 0);
        return freePoleInfoMap[p];
    }

    bool _in_set( set<VertexID > &s, VertexID v ){
        return ( s.count(v) > 0 );
    }
    
    void MainStructure::reAlignIDs( HMesh::VertexIDRemap &remapper ){
        
        PoleInfoMap p;
        freePolesSet.clear();
        
        for( int i = 0; i < freePoles.size(); ++i ){
            VertexID newID = remapper[freePoles[i]];
            p[newID] = freePoleInfoMap[freePoles[i]];
            freePoles[i] = newID;
            freePolesSet.insert( newID );
        }
        for( int i = 0; i < gluedPoles.size(); ++i ){
            gluedPoles[i] = remapper[gluedPoles[i]];
        }
        freePoleInfoMap = std::move( p );
        skel->reAlignIDs( remapper );
    }
    
    void MainStructure::glueModule( Module &m, vector<Match> &matches ){
        set<VertexID> glued_m_poles;
        set<VertexID> glued_h_poles;
        
        // for assert purposes
        size_t  old_free_poles_size = freePoles.size(),
                old_glued_poles_size = gluedPoles.size();
        
        for( Match match : matches){
            assert( find( m.poleList.begin(), m.poleList.end(), match.first ) != m.poleList.end( ));
            assert( find( freePoles.begin(), freePoles.end(), match.second ) != freePoles.end( ));
            glued_m_poles.insert( match.first );
            glued_h_poles.insert( match.second );
        }
        
        for( VertexID v : m.poleList ){
            if( glued_m_poles.count(v) > 0 ){
                gluedPoles.push_back( v );
            }
            else{
                freePoles.push_back( v );
                freePolesSet.insert( v );
                freePoleInfoMap[v] = m.getPoleInfo( v );
            }
        }
        // remove from freePoles the host poles involved  and put them into gluedPoles
        for( VertexID v : glued_h_poles ){
            freePoles.erase( remove(freePoles.begin(), freePoles.end(), v));
            freePolesSet.erase( v );
            gluedPoles.push_back( v );
            assert( freePoleInfoMap.count(v) > 0 );
            freePoleInfoMap.erase( v );
        }
        assert( glued_h_poles.size() == glued_m_poles.size() );
        assert( freePoles.size() == freePolesSet.size());
        assert( old_glued_poles_size + glued_m_poles.size() + glued_h_poles.size() == gluedPoles.size() );
        assert( old_free_poles_size + ( m.poleList.size() - glued_m_poles.size( ) - glued_h_poles.size())
                == freePoles.size());
        assert( old_glued_poles_size + old_free_poles_size + m.poleList.size()
                == freePoles.size() + gluedPoles.size());
        
        GluedModuleInfo gmi;
        gmi.module              = &m;
        gmi.t_start             = time;
        gmi.connection_valence  = matches.size();

        modules.push_back( gmi );
        ++time;
        
        /***** DEBUG AND SANITY CHECK ****/
        cout << glued_m_poles.size() << "-valent glueing at time : " << time << endl;
        cout << " num of free poles " << freePoles.size() << " # set : " << freePolesSet.size();
        cout << " num of glued poles " << gluedPoles.size() << endl;
        /*****          END         ****/
        
        matches.clear();
        skel->merge( m.getSkeleton(), matches );
        skel->saveToFile( "//Users//francescousai//Desktop//example.skel" );
        skel->saveCollisionDetectionHierarchyToFile( "//Users//francescousai//Desktop//example_cd.skel" );
    }
    
    
    
    
    bool MainStructure::isColliding(const Module &m) const{
        bool collision_found = false;
        
#ifdef TRACE
        cout << "module sphere : " << m.bsphere_center << " # " << m.bsphere_radius << endl;
#endif
        
        for ( size_t i = 0; i < modules.size() && (!collision_found); ++i) {
            Module *mi = modules[i].module;
            cout << "gluedM sphere : " << mi->bsphere_center << " # " << mi->bsphere_radius << endl;
            
            double collision_deepness =
                sphere_intersects( m.bsphere_center, m.bsphere_radius,
                                   mi->bsphere_center, mi->bsphere_radius);
            double min_radius = min( mi->bsphere_radius, m.bsphere_radius );
            collision_found = collision_deepness > ( min_radius / 2.0 );
            if( collision_found ){
                cout << "collision found! " << collision_deepness << " => " << min_radius << endl;
            }

        }
        return collision_found;
    }
    
}