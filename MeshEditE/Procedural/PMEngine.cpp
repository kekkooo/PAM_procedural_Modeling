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
#include <set>

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
        if(m)
        {
            CGLA::gel_srand( m->no_vertices()) ;
            CGLA::gel_srand( m->no_faces() );
            CGLA::gel_srand( m->no_halfedges() );
        }
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
//        Procedural::Operations::create_basic_PAM( *m, 0.5 );
        Procedural::Operations::create_PAM_box( *m, 3.0, 1.0, 2.0 );
        _polesList.Update( m );
    }
    
    void Engine::setMesh(HMesh::Manifold *mesh)
    {
        m = mesh;
        invalidateAll();
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
        
        CGLA::Vec3d axes[6];
        axes[0] = CGLA::Vec3d(  0.0,  1.0,  0.0 );
        axes[1] = CGLA::Vec3d(  0.0, -1.0,  0.0 );
        axes[2] = CGLA::Vec3d(  1.0,  0.0,  0.0 );
        axes[3] = CGLA::Vec3d( -1.0,  0.0,  0.0 );
        axes[4] = CGLA::Vec3d(  0.0,  0.0,  1.0 );
        axes[5] = CGLA::Vec3d(  0.0,  0.0, -1.0 );
        
        for( VertexID pole_id : poles)
        {
            if( trajectories.count( pole_id ) > 0 )
            {
                bool shake = trajectories[pole_id].total_length > _geometric_info.MeanLength() * _polesList.PoleAge( pole_id );
                trajectories[pole_id].no_calls = ( trajectories[pole_id].no_calls + 1 );// % ( CGLA::gel_rand() % 100 + 10 );
                
                if( trajectories[pole_id].no_calls == 1 || shake )
                {
                    trajectories[pole_id].d1 = CGLA::gel_rand() % 6;
                    trajectories[pole_id].d2 = ((CGLA::gel_rand() % 6) + 4 ) % 6;
                    trajectories[pole_id].current_dir = vertex_normal( *m, pole_id );
                    trajectories[pole_id].current_dir.normalize();
                    trajectories[pole_id].total_length = 0.0;
                }
                
                CGLA::Vec3d dir = ( trajectories[pole_id].current_dir
                                  + axes[trajectories[pole_id].d1] * cos( 0.05 * (double)trajectories[pole_id].no_calls )
                                  + axes[trajectories[pole_id].d2] * sin( 0.05 * (double)trajectories[pole_id].no_calls ));
                dir.normalize();
                dir *= _geometric_info.MeanLength();
                
                trajectories[pole_id].total_length += dir.length();
                
                Geometric::extrude_pole(*m, pole_id, dir, true, 0.9);
                trajectories[pole_id].no_calls++;
            }
            else
            {
                Geometric::extrude_pole( *m, pole_id, NAN, true, 0.9 );

                trajectories[pole_id].current_dir = vertex_normal( *m, pole_id );
                trajectories[pole_id].current_dir.normalize();
                trajectories[pole_id].total_length = _geometric_info.MeanLength();
                trajectories[pole_id].d1 = 0;
                trajectories[pole_id].d1 = 1;
                trajectories[pole_id].no_calls = 1;
            }
        }
        m->cleanup();
        _edges_info_container.Update( m, true, true );

    }
    
    void Engine::polarSubdivision( )
    {
        polar_subdivide( *m, 1 );
        m->cleanup();
        _edges_info_container.Update( m, true, true );

    }
    
    void Engine::perturbate()
    {
        _polesList.Update( m, true);
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update( m, _edges_info_container );
        
        vector< VertexID >  selected;
        vector< double >    per_vertex_ratio;
//        int                 distance = (int) ( log2( _polesList.No_Poles( )) + log2( _polesList.MeanPoleValency( )));
        int                 distance = (int)log2( _polesList.MeanPoleValency( ));
        // ratio should be calculeted per vertex, accordingly to the region rib_edge_loop radius
        double              ratio    = 0.2;
        double              avg_length   = _geometric_info.MeanLength();
        cout << "adding noise with distance " << distance << endl;
        
        for( VertexID vid : m->vertices() )
        {
            assert(vid != InvalidVertexID);

            int vertex_distance = _geometric_info.CombinedDistance()[vid];
            if ( vertex_distance > distance )
            {
                double radius;
                map< VertexID, double > radii;
                // here is possible to optimize things!
                double mean_radius = ring_mean_radius( *m, get_first_rib_edge
                                                      ( *m, vid, _edges_info_container.edgeInfo( )),
                                                      radii );
                radius = radii[vid];
                double alt_ratio   = ( radius * (log2( (double)vertex_distance )) * _geometric_info.MeanLength() )  ;
//                double alt_ratio = _geometric_info.MeanLength() / mean_radius;
//                cout << " alt ratio is : " << alt_ratio << " avg_radius : " << mean_radius << " edge_avg : " << _geometric_info.MeanLength()
//                     << "dist : " << vertex_distance << " inverse_log_dist " << ( 1.0 / log2(vertex_distance)) << endl;
                selected.push_back( vid );
                per_vertex_ratio.push_back( alt_ratio );
            }
        }
//        add_perpendicular_noise( *m, selected, ratio, avg_length );
        add_perpendicular_noise( *m, selected, per_vertex_ratio, avg_length );
    }
   
    
    void Engine::smooth_near_junctions()
    {
        _polesList.Update( m, true);        
        if (_polesList.No_Poles() <= 2 ) return; // there are no junctions
        
        
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update( m, _edges_info_container );
        
        vector< VertexID >  selected;
        //        size_t              distance = ( m->no_vertices() ) / ( _polesList.MeanPoleValency() * _polesList.No_Poles());
        int              distance = (int) ( log2( _polesList.No_Poles( )) + log2( _polesList.MeanPoleValency( )));
        // ratio should be calculeted per vertex, accordingly to the region rib_edge_loop radius
        cout << "adding noise with distance " << distance << endl;
        
        for( VertexID vid : m->vertices() )
        {
            assert(vid != InvalidVertexID);
            if ( _geometric_info.CombinedDistance()[vid] < distance && !_polesList.IsPole( vid ))
            {
                selected.push_back( vid );
            }
        }
        //        add_noise( *m, VertexType::REGULAR, ratio, cutoff, selected );
        Procedural::Operations::Algorithms::selected_vertices_inverse_distance_laplacian(*m, selected);

    }
    
    
    
// HERE ARE DEFINED SOME FUNCTIONS THAT MUST BE DELETED WHEN THE THING IS WORKING PROPERLY
    void Engine::addRandomBranches()
    {
        _polesList.Update(m);
        _edges_info_container.Update( m, true, true );
        _geometric_info.Update(m, _edges_info_container, true );
        
        int max_new_branches = (int) log2(_polesList.No_Poles( ));
        
        CGLA::gel_rand();
        CGLA::gel_rand();
        CGLA::gel_rand();
        
        std::cout << "there are  " << _polesList.No_Poles() << " branches " << endl;
        // save current poles, in order to find which were added.
        size_t  old_size = _polesList.No_Poles();
        bool there_are_junctions = _polesList.No_Poles() > 2;

        int count = 0;
        for( VertexID vid : m->vertices() )
        {
            bool branch_added = false;
            if( _polesList.IsPole( vid )) continue;
            int p = CGLA::gel_rand() % 100;
            
            cout << " p is : " << p << endl;
            
            if( p == 50 ) // it does not mean too much but should work
            {
                CGLA::Vec3d up( 0.0, 1.0, 0.0 );
                // add branches only on the upper side
                double angle = Geometry::angle( up, vertex_normal(*m, vid) );
                if( angle > M_PI_2 ) continue;

                
                
                int distance_limit   = (int) ( log2( _polesList.No_Poles( )) + log2( _polesList.MeanPoleValency( )));
                int branch_size      =  ( CGLA::gel_rand() % 3 ) + 2; // this should depend on the branch thickness
                
                cout << _polesList.No_Poles() << " branches. mean valence is : " << _polesList.MeanPoleValency() <<
                     " you are using distance limit : " << distance_limit << " with branch size : " << branch_size << endl;
                
//                if( _geometric_info.CombinedDistance()[vid].first > distance_limit )
                if( there_are_junctions )
                {
                    
                    if( _geometric_info.JunctionDistance()[vid] > distance_limit
                     && _geometric_info.PoleDistance()[vid] > branch_size * 2 )
                    {
                        buildCleanSelection();
                        add_branch( *m, vid, branch_size, vertex_selection );
                        ++count;
                        branch_added = true;
                    }
                }
                else
                {
                    if( _geometric_info.PoleDistance()[vid] > distance_limit * branch_size )
                    {
                        buildCleanSelection();
                        add_branch( *m, vid, branch_size, vertex_selection );
                        ++count;
                        branch_added = true;
                    }
                }
                if( branch_added )
                {
                    _polesList.Update(m);
                    _edges_info_container.Update( m, true, true );
                    _geometric_info.Update(m, _edges_info_container, true );
                }
            }
            if ( count >= max_new_branches ) break;
        }
        

        if( count > 0)
        {
            std::cout << "added " << count << " branches " << endl;
            _polesList.Update( m, true );
            std::cout << "now there are  " << _polesList.No_Poles() << " branches " << endl;
            
            assert( _polesList.No_Poles() == old_size + count );
            
            for( VertexID pole : _polesList.Poles() )
            {
                // flatten only poles with age 0
                if( _polesList.PoleAge( pole ) == 0)
                    flatten_pole( *m, pole );
            }
            
            m->cleanup();
        }
    }
    
    // should have be done better,
    // but I need a proportional scaling along the branches' rib edge loops
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