//
//  PMEngine.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 30/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "PMEngine.h"
#include <MeshEditE/Procedural/Operations/basic_shapes.h>
#include <MeshEditE/Procedural/Operations/Algorithms.h>
#include <MeshEditE/Procedural/Operations/geometric_operations.h>

using namespace Procedural::Operations;

namespace Procedural
{
    Engine::Engine()
    {
    }

    void Engine::buildCube(Manifold &mesh)
    {
        invalidateAll();
        Procedural::Operations::create_basic_PAM( mesh, 0.5 );
        _polesList.Update( m );
    }

    void Engine::buildCube()
    {        
        invalidateAll();
        Procedural::Operations::create_basic_PAM( *m, 0.5 );
        _polesList.Update( m );
    }
    
    void Engine::setMesh(HMesh::Manifold *mesh)
    {
        m = mesh;
        _polesList.Update( m, true );
        _edges_info_container.Update( m );
    }
    
    
    void Engine::smooth( Manifold& mesh )
    {
        Procedural::Operations::Algorithms::cotangent_weights_laplacian_smoothing( mesh );
    }
    
    void Engine::extrudePoles()
    {
        // REMEMBER THAT MODIFICATIONS OF THE MESH OUTSIDE THE ENGINE WILL NOT
        // REVEAL CHANGES ON polesLIST nor edges_info_container so you need to update
        if( !m ) return;
//        assert( _polesList.IsValid() );

        //need to remove this later
        _polesList.Update(m);
        
        auto poles = _polesList.Poles();
        invalidateAll();
        
        for( VertexID vid : poles)
        {
            Geometric::extrude_pole( *m, vid, NAN, true, 0.8 );
        }

        _edges_info_container.Update( m, true, true );

        if( poles.size() > 2 )
        {
            for( VertexID vid : poles )
            {
                // here there is a problem because if there are only 2 poles
                // there will be no junctions and smooth_pole "walks"
                // until it finds a junction. Hence if poles.size() == 2 => infinite loop
                Geometric::smooth_pole( *m, vid, _edges_info_container.edgeInfo( ));
            }
        }
//        _polesList.Update( m );
    }
    
    
    void Engine::polarSubdivision( )
    {
        polar_subdivide( *m, 1 );
    }
    
}