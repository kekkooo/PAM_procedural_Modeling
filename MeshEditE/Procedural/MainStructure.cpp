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

    const PoleList& MainStructure::getPoleList(){
        return freePoles;
    }
    
    bool _in_set( set<VertexID > &s, VertexID v ){
        return (s.count(v) > 0);
    }
    
    void MainStructure::glueModule( Module &m, vector<Match> &matches ){
        set<VertexID> glued_m_poles;
        set<VertexID> glued_h_poles;
        
        // for assert purposes
        int old_free_poles_size = freePoles.size(),
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
            }
        }
        // remove from freePoles the host poles involved  and put them into gluedPoles
        for( VertexID v : glued_h_poles ){
            freePoles.erase( remove(freePoles.begin(), freePoles.end(), v));
            gluedPoles.push_back(v);
        }
        
        assert( old_glued_poles_size + glued_m_poles.size() == gluedPoles.size() );
        assert( old_free_poles_size + ( m.poleList.size() - glued_m_poles.size( )) == freePoles.size());
        assert( old_glued_poles_size + old_free_poles_size + m.poleList.size()
                == freePoles.size() + gluedPoles.size());
        
        GluedModuleInfo gmi;
        gmi.module              = &m;
        gmi.t_start             = time;
        gmi.connection_valence  = matches.size();

        modules.push_back( gmi );
    }
}