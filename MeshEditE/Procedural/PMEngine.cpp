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
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <GEL/CGLA/CGLA.h>
#include <MeshEditE/Test.h>

using namespace Procedural::Operations;
using namespace Procedural::Operations::Geometric;
using namespace Procedural::Operations::Structural;
using namespace Procedural::Structure;
using namespace Procedural::Geometry;

namespace Procedural
{
    Engine::Engine()
    {
        invalidateAll();
    }
    
    void Engine::invalidateAll()
    {
        invalidateEdgeInfo();
        invalidatePolesList();
        invalidateGeometricInfo();
    }
    
    void Engine::buildCleanSelection()
    {
        vertex_selection.clear();
        for( auto v : m->vertices( )) { vertex_selection[v] = 0; }
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
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update( m, _edges_info_container );
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
        _polesList.Update( m, true );
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update( m, _edges_info_container );
        
        auto poles = _polesList.Poles();
        invalidateAll();
        
        CGLA::Vec3d up( 0.0, 1.0, 0.0), down( 0.0, -1.0, 0.0 );
        
        for( VertexID pole_id : poles)
        {
            if( trajectories.count( pole_id ) > 0 )
            {
                bool au_contraire = trajectories[pole_id].no_calls < 15;
                trajectories[pole_id].no_calls = ( trajectories[pole_id].no_calls + 1 );// % ( CGLA::gel_rand() % 100 + 10 );
                
                CGLA::Vec3d dir;
//                CGLA::Vec3d v_normal = vertex_normal( *m, pole_id );
                double angle = Geometry::angle( up, trajectories[pole_id].current_dir );
                if( angle < M_PI_2 || ( angle > M_PI_2 && au_contraire) )
                {
//                    dir = ( v_normal + up * log(trajectories[pole_id].no_calls++ ) );
                    dir = ( trajectories[pole_id].current_dir + up * cos( 0.05 * (double)trajectories[pole_id].no_calls++ ) );
                    dir.normalize();
                    dir *= _geometric_info.MeanLength();
                }
                else
                {
//                    dir = ( v_normal + down * log(trajectories[pole_id].no_calls++ ) );
                    dir = ( trajectories[pole_id].current_dir + down * sin( 0.05 * (double) trajectories[pole_id].no_calls++ ));
                    dir.normalize();
                    dir *= _geometric_info.MeanLength();
                }
                
                trajectories[pole_id].total_length += dir.length();
                
                Geometric::extrude_pole(*m, pole_id, dir, true, 0.9);
            }
            else
            {
                Geometric::extrude_pole( *m, pole_id, NAN, true, 0.9 );

                trajectories[pole_id].current_dir = vertex_normal( *m, pole_id );
                trajectories[pole_id].current_dir.normalize();
                trajectories[pole_id].total_length = _geometric_info.MeanLength();
                trajectories[pole_id].no_calls = 1;
            }
        }

        _edges_info_container.Update( m, true, true );
    }
    
    void Engine::polarSubdivision( )
    {
        polar_subdivide( *m, 1 );
        _edges_info_container.Update( m, true, true );
    }
    
    void Engine::perturbate()
    {
        _polesList.Update(m);
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update( m, _edges_info_container );
        
        vector< VertexID >  selected;
        int                 distance = m->no_vertices() / ( _polesList.MeanPoleValency() * _polesList.No_Poles());
        double              ratio    = 0.2;
        double              cutoff   = _geometric_info.MeanLength();
        cout << "adding noise with distance " << distance << endl;
        
        for( VertexID vid : m->vertices() )
        {
            if ( _geometric_info.CombinedDistance()[vid].first > distance )
            {
                selected.push_back( vid );
            }
        }
        add_noise( *m, VertexType::REGULAR, ratio, cutoff, selected );
    }
    
    
    
// HERE ARE DEFINED SOME FUNCTIONS THAT MUST BE DELETED WHEN THE THING IS WORKING PROPERLY
    void Engine::addRandomBranches()
    {
        _polesList.Update(m);
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update(m, _edges_info_container, true );
        
        int distance_limit  = 10;
        int branch_size     =  ( CGLA::gel_rand() % 3 ) + 2;
        
        CGLA::gel_srand(0);
        CGLA::gel_rand();
        CGLA::gel_rand();
        CGLA::gel_rand();
        
        std::cout << "there are  " << _polesList.No_Poles() << " branches " << endl;
        int count = 0;
        for( VertexID vid : m->vertices() )
        {
            if( _polesList.IsPole( vid )) continue;
            int p = CGLA::gel_rand() % 200;
            if( p == 100 ) // it does not mean too much but should work
            {
//                if( _geometric_info.CombinedDistance()[vid].first > distance_limit )
                if( _geometric_info.CombinedDistance()[vid].first > distance_limit )
                {
                    buildCleanSelection();
                    add_branch( *m, vid, branch_size, vertex_selection );
                    ++count;
                }
            }
            if ( count >= 5 ) break;
        }
        
        std::cout << "added " << count << " branches " << endl;
        _polesList.Update( m, true );
        for( VertexID pole : _polesList.Poles() )
            flatten_pole( *m, pole );
        
        std::cout << "now there are  " << _polesList.No_Poles() << " branches " << endl;
        
    }
    void Engine::pickABranchAndScaleIt( int mode )
    {
        _polesList.Update( m );
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update( m, _edges_info_container );
        
        for( auto pole : _polesList.Poles())
        {
            
//            VertexID pole =
//                _polesList.Poles()[ CGLA::gel_rand(mode * CGLA::gel_rand())
//                                    % _polesList.No_Poles() ];
            
            vector<HalfEdgeID> ring_starters;
            bool there_are_junctions = _polesList.No_Poles() > 2;
            
            if( there_are_junctions )
            {
                get_rings_from_pole( *m, pole,_edges_info_container.edgeInfo(), ring_starters );
            }
            else
            {
                //calculate a reasonable max_iterations
                int max_iter = m->no_vertices() / ( valency( *m, pole ) * _polesList.No_Poles() );
                get_rings_from_pole( *m, pole, _edges_info_container.edgeInfo(), ring_starters, max_iter );
            }
            
            cout << "scaling in mode " << mode << endl;
            
            switch (mode) {
                case 0:
                {
                    for( auto ring : ring_starters )
                    {   scale_ring_radius( *m, ring, 1.05, true );    }
                }
                    break;
                    
                case 1 :
                {
                    for( auto ring : ring_starters )
                    {    scale_ring_radius( *m, ring, 0.95, true );   }
                }
                    break;
                case 2 :
                {
                    vector< double > scale_factors;
                    linspace( 1.01, 1.15, ring_starters.size(), scale_factors );

                    assert(scale_factors.size() == ring_starters.size());
                    for( int i = 0; i < ring_starters.size(); ++i )
                    {
                        scale_ring_radius( *m, ring_starters[i], scale_factors[i], true );
                    }
                }
                case 3 :
                {
                    vector< double > scale_factors;
                    linspace( 1.15, 1.01, ring_starters.size(), scale_factors );
                    
                    assert(scale_factors.size() == ring_starters.size());
                    for( int i = 0; i < ring_starters.size(); ++i )
                    {
                        scale_ring_radius( *m, ring_starters[i], scale_factors[i], true );
                    }
                }

                    break;
                default:
                    break;
            }
        }

        

    }
    
    
// STOPS HERE
}