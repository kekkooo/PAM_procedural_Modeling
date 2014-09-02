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
        
        for( VertexID pole_id : poles)
        {
            // find matching vertex # naive for now
            VertexID matching = poles[ CGLA::gel_rand(pole_id.get_index()) % poles.size() ];

            // find a matching vertex that is far from it
//            VertexID matching = _pole_tracking.GetMatch( pole_id );
//                
//            if( matching == InvalidVertexID )
//            {
//                int distance = CGLA::gel_rand() % m->no_vertices() / 2;
//                Walker w = m->walker( pole_id );
//                while( distance-- > 0 )
//                {
//                    w = w.next().opp().next();
//                }
//                matching = w.vertex();
//                _pole_tracking.AddMatch( pole_id, matching );
//            }
            assert( matching != InvalidVertexID );

            if( matching != pole_id )
            {
                CGLA::Vec3d dir =
                    ( m->pos( pole_id ) - m->pos( matching ) )
                +   ( Geometry::vertex_normal(*m, pole_id ) * (double)( CGLA::gel_rand() % 5 ) );
//                +   ( Geometry::vertex_normal(*m, matching_pole));
                dir.normalize();
                double dir_mult = std::max( 0.5 * _geometric_info.MeanLength(), mean_length_of_outoing_he( *m, pole_id));
                dir *= dir_mult;
                Geometric::extrude_pole(*m, pole_id, dir, true, 0.8);
            }
            else
            {
                Geometric::extrude_pole( *m, pole_id, NAN, true, 0.8 );
            }
        }

        _edges_info_container.Update( m, true, true );

//        if( poles.size() > 2 )
//        {
//            for( VertexID vid : poles )
//            {
//                // here there is a problem because if there are only 2 poles
//                // there will be no junctions and smooth_pole "walks"
//                // until it finds a junction. Hence if poles.size() == 2 => infinite loop
//                Geometric::smooth_pole( *m, vid, _edges_info_container.edgeInfo( ));
//            }
//        }
//        _polesList.Update( m );
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
        int                 distance = m->no_vertices() / ( _polesList.MeanPoleValency() * 8);
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
        _geometric_info.Update(m, _edges_info_container );
        
        int distance_limit  = 10;
        int branch_size     = 3;
        
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
                if( _geometric_info.PoleDistance()[vid].first > distance_limit )
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
        
        VertexID pole =
            _polesList.Poles()[ CGLA::gel_rand(mode * CGLA::gel_rand())
                                % _polesList.No_Poles() ];
        
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
                {   scale_ring_radius( *m, ring, 1.05 );    }
            }
                break;
                
            case 1 :
            {
                for( auto ring : ring_starters )
                {    scale_ring_radius( *m, ring, 0.95 );   }
            }
                break;
            case 2 :
            {
                vector< double > scale_factors;
                linspace( 1.01, 1.15, ring_starters.size(), scale_factors );

                assert(scale_factors.size() == ring_starters.size());
                for( int i = 0; i < ring_starters.size(); ++i )
                {
                    scale_ring_radius( *m, ring_starters[i], scale_factors[i] );
                }
            }
            case 3 :
            {
                vector< double > scale_factors;
                linspace( 1.15, 1.01, ring_starters.size(), scale_factors );
                
                assert(scale_factors.size() == ring_starters.size());
                for( int i = 0; i < ring_starters.size(); ++i )
                {
                    scale_ring_radius( *m, ring_starters[i], scale_factors[i] );
                }
            }

                break;
            default:
                break;
        }

        

    }
    
    
// STOPS HERE
}