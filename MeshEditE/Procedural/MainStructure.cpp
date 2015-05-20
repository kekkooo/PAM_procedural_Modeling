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

using namespace std;
using namespace HMesh;
using namespace Procedural;
using namespace Procedural::GraphMatch;

namespace Procedural{
    
    MainStructure::MainStructure(){
        time = 0;
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

    bool _in_set( set<VertexID > &s, VertexID v ){
        return (s.count(v) > 0);
    }
    
    void MainStructure::reAlignIDs( HMesh::VertexIDRemap &remapper ){
        
        freePolesSet.clear();
        for( int i = 0; i < freePoles.size(); ++i ){
            VertexID newID = remapper[freePoles[i]];
            freePoles[i] = newID;
            freePolesSet.insert( newID );
        }
        for( int i = 0; i < gluedPoles.size(); ++i ){
            gluedPoles[i] = remapper[gluedPoles[i]];
        }
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
            }
        }
        // remove from freePoles the host poles involved  and put them into gluedPoles
        for( VertexID v : glued_h_poles ){
            freePoles.erase( remove(freePoles.begin(), freePoles.end(), v));
            freePolesSet.erase( v );
            gluedPoles.push_back(v);

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


        

    }
}