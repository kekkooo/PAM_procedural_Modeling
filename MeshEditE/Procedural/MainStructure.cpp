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

#include<MeshEditE/Procedural/Helpers/geometric_properties.h>

#include "Test.h"

using namespace std;
using namespace HMesh;
using namespace Procedural;
using namespace Procedural::Geometry;

namespace Procedural{
    
    MainStructure::MainStructure(){
        time = 0;
        skel = new Skeleton();
    }

    const PoleList& MainStructure::getPoles() const {
        return freePoles;
    }

    const PoleList& MainStructure::getFreePoles() const{
        return freePoles;
    }
    
    const PoleSet& MainStructure::getFreePoleSet() const{
        return freePolesSet;
    }

    const PoleList& MainStructure::getGluedPoles() const{
        return freePoles;
    }
    
    const PoleInfo& MainStructure::getPoleInfo( HMesh::VertexID p ) const{
        assert(freePoleInfoMap.count(p) > 0);
        return freePoleInfoMap.at( p );
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
    
    void MainStructure::glueModule( const HMesh::Manifold& mani, Module &m, vector<Match> &matches ){
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
                freePoleInfoMap[v]              = m.getPoleInfo( v );
                freePoleInfoMap[v].geometry.pos = mani.pos( v );
                CGLA::Vec3d normal = vertex_normal( mani, v );
                normal.normalize();
                freePoleInfoMap[v].geometry.normal = normal ;
                

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
        /*
        for( const Match& match : matches ){
            cout << "( " << match.first << ", " << match.second << " )  # ";
        }
         */
        cout << endl;
        cout << " num of free poles " << freePoles.size() << " # set : " << freePolesSet.size();
        cout << " num of glued poles " << gluedPoles.size() << endl;
        /*****          END         ****/
        
        matches.clear();

        skel->merge( m.getSkeleton(), matches );
//        skel->saveToFile( "//Users//francescousai//Desktop//example.skel" );
        

    }
    
    bool MainStructure::isColliding(const Module &m) const{
        return Procedural::collide( *( this->skel ), m.getSkeleton() );
    }
    
    void MainStructure::saveBVH( std::string path ) const{
        skel->saveCollisionDetectionHierarchyToFile( path );
    }
    
    void MainStructure::saveSkeleton( std::string path ) const{
        skel->saveToFile( path );
    }
    
    void MainStructure::setBoundingSphere( CGLA::Vec3d center, double radius ){
        skel->bounding_sphere.center = center;
        skel->bounding_sphere.radius = radius;
    }
    
    
}