//
//  StatefulEngine.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 17/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "StatefulEngine.h"

#include "polarize.h"

#include "MeshEditE/Procedural/Helpers/manifold_copy.h"
#include "MeshEditE/Procedural/Helpers/geometric_properties.h"


using namespace std;

using namespace HMesh;
using namespace CGLA;

using namespace Procedural::Engines;
using namespace Procedural::Helpers;
using namespace Procedural::Geometry;
using namespace Procedural::Helpers::ModuleAlignment;
using namespace Procedural::GraphMatch;

/*=========================================================================*
 *                     PRIVATE FUNCTIONS                                   *
 *=========================================================================*/

StatefulEngine::StatefulEngine()
{
    this->m = NULL;
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    randomizer.seed( seed );
    treeIsValid = false;
    dim_constraint = DimensionalityConstraint::Constrained_3D;
}


void StatefulEngine::buildRandomTransform( CGLA::Mat4x4d &t ){
    // must be initialized
    float rand_max =  static_cast<float>( randomizer.max( )); //RAND_MAX / ( M_PI * 2.0 );
    float x1 = static_cast<float>( randomizer( )) / rand_max,
    x2 = static_cast<float>( randomizer( )) / rand_max,
    x3 = static_cast<float>( randomizer( )) / rand_max;
    
    t = ModuleAlignment::random_rotation_matrix_arvo( x1, x2, x3 );
}


void StatefulEngine::buildCollisionAvoidingTranslation( const CGLA::Mat4x4d &rot, CGLA::Mat4x4d &tr )
{
    vector< Vec3d > transformed_points;
    Vec3d  H_centroid, M_centroid;
    double H_radius, M_radius;
    // save transformed module vertices
    for( VertexID vid : M_vertices ){
        transformed_points.push_back( rot.mul_3D_point( this->m->pos( vid )));
    }
    // find its bsphere and centroid
    bsphere( transformed_points, M_centroid, M_radius );
    bsphere( *(this->m), H_vertices, H_centroid, H_radius );
    // build the collision avoiding translation
    Vec3d   dir     = M_centroid - H_centroid;
    double  length  = H_radius + M_radius - dir.length();
    dir.normalize();
    dir *= length;
    
    tr = CGLA::translation_Mat4x4d( dir );
}


void StatefulEngine::buildCollisionAvoidingRandomTransform( CGLA::Mat4x4d &t ){
    Mat4x4d rot, tr;
    // build random transform
    buildRandomTransform( rot );
    // build collision avoiding translation
    buildCollisionAvoidingTranslation( rot, tr );
    // return composition
    t = tr * rot;
}


void StatefulEngine::buildHostKdTree(){
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    
//    ModuleAlignment::build_manifold_kdtree( (*this->m), this->H_vertices, this->tree );
    ModuleAlignment::build_manifold_kdtree( (*this->m), this->H_candidates.getCandidates(), this->tree );
    treeIsValid = true;
}


void StatefulEngine::transformModulePoles( CGLA::Mat4x4d &t, VertexPosMap &new_pos ){
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    assert( this->M_vertices.size() > 0 );
    
    for( VertexID vid : M_vertices ){
        if( !is_pole( (*this->m), vid )) continue;
        new_pos[vid] = t.mul_3D_point( m->pos( vid ));
    }
}


bool ID_and_dist_comparer( ID_and_dist &l, ID_and_dist &r ) { return ( l.second < r.second ); }


void StatefulEngine::matchModuleToHost( VertexPosMap& module_poles_positions, VertexMatch& M_pole_to_H_vertex ){
    set<VertexID> assigned_vertices;
    for( auto id_and_pos : module_poles_positions)
    {
        assert( is_pole( *this->m, id_and_pos.first ));
        double      distance    = numeric_limits<double>::max();
        Vec3d       foundVec;
        VertexID    foundID = InvalidVertexID, second_choiceID = InvalidVertexID;
        // here there could be problems because I should be sure that there will not be
        // two poles assigned to the same candidate vertex on the host
        // closest point should always return true
        bool        have_found  = tree.closest_point( id_and_pos.second , distance, foundVec, foundID );
        
        cout << id_and_pos.first << " # " << id_and_pos.second << "# matched : " << foundID << " # dist : " << distance << endl;
        assert( have_found );
        assert( foundID != InvalidVertexID );
        
        if( assigned_vertices.count( foundID ) == 0 )
        {
            assigned_vertices.insert( foundID );
            M_pole_to_H_vertex[id_and_pos.first] = foundID;
        }
        else
        {
            findSecondClosest( foundID, second_choiceID, assigned_vertices );
            assert( second_choiceID != InvalidVertexID );
            assigned_vertices.insert( second_choiceID );
            M_pole_to_H_vertex[id_and_pos.first] = second_choiceID;
        }
    }
}


void StatefulEngine::findSecondClosest(const HMesh::VertexID &closest, HMesh::VertexID &second_closest, VertexSet &assigned){
    assert( assigned.count(closest) > 0 );
    // consider closest
    Vec3d       cp          = m->pos( closest );
    // find the nearest of the vertices in its 1-ring
    double      min_dist    = numeric_limits<double>::max();
    for( Walker w = m->walker( closest ); !w.full_circle(); w = w.circulate_vertex_ccw())
    {
        if( assigned.count(w.vertex())) { continue; }   // skip all assigned vertices
        double dist = ( cp - m->pos( w.vertex( ))).length( );
        if( dist < min_dist ){ min_dist = dist; }
    }
    // take that distance as radius and use it with tree.in_sphere with center in closest
    double              radius          = min_dist * 1.1;
    vector<VertexID>    in_sphere_ID;
    vector<Vec3d>       in_sphere_points;
    IDs_and_dists       ids_and_dists;
    
    tree.in_sphere( cp, radius, in_sphere_points, in_sphere_ID );
    // take the nearest point
    assert( in_sphere_ID.size() == in_sphere_points.size());
    for( int i = 0; i < in_sphere_points.size(); ++i)
    {
        ids_and_dists.push_back( make_pair(in_sphere_ID[i], (cp - in_sphere_points[i]).length()));
    }
    // sort in ascending order accordingly to the distance
    std::sort( ids_and_dists.begin(), ids_and_dists.end(), ID_and_dist_comparer );
    IDs_and_dists::iterator it      = ids_and_dists.begin();
    bool                    done    = ( assigned.count( it->first ) == 0 );
    while ( !done )
    {
        ++it;
        assert( it != ids_and_dists.end( ));
        done    = ( assigned.count( it->first ) == 0 );
    }
    second_closest = it->first;
    assigned.insert( it->first );
}


void StatefulEngine::alignUsingBestMatch( ){
    if( best_match.IsValid( )){
        apply_optimal_alignment( *m, M_vertices, best_match.getMatchInfo() );
        align_module_normals_to_host( *m, M_vertices, best_match.getMatchInfo().matches );
    }
}

void StatefulEngine::actualGlueing(){
    if( best_match.IsValid( )){
        add_necessary_poles( *m, best_match.getMatchInfo().matches, H_vertices, M_vertices );
        glue_matches( *m, best_match.getMatchInfo().matches );
        consolidate();
    }
}

void StatefulEngine::fillCandidateSet(){
#warning nowadays it considers just the poles
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    
    for( VertexID vid : m->vertices( )){
        if( is_pole( *m, vid )) {
            CandidateInfo _;
            H_candidates.insert( vid, _ );
        }
    }
}

/*=========================================================================*
 *                     PUBLIC FUNCTIONS                                    *
 *=========================================================================*/

StatefulEngine& StatefulEngine::getCurrentEngine(){
    static StatefulEngine instance;
    return instance;
}


void StatefulEngine::setHost( Manifold &host ){
//    assert( this->m == NULL );
    this->m = &host;
}

void StatefulEngine::setModule( Manifold &module ){
    assert( this->m != NULL );
    assert( this->H_vertices.size() == 0 );
    assert( this->M_vertices.size() == 0 );
    add_manifold( (*this->m), module, H_vertices, M_vertices );
    fillCandidateSet();
    buildHostKdTree();    
}

void StatefulEngine::consolidate(){
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    assert( this->M_vertices.size() > 0 );
    // in this way you lose any reference to which vertices are from host and module
    H_vertices.clear();
    M_vertices.clear();
    H_candidates.clear();
    treeIsValid = false;
}


void StatefulEngine::testMultipleTransformations( int no_tests, int no_glueings ){
        assert( this->m != NULL );
        assert( this->H_vertices.size() > 0 );
        assert( this->M_vertices.size() > 0 );
    
    assert( no_tests > 0 );
    
    vector< matches_and_cost >    matches_vector;
    vector< match_info >          proposed_matches;
    
    for( int i = 0; i < no_tests; ++i ){
        VertexPosMap    transformed_M_poles;
        vertex_match    M_to_H;
        vector<Match>   current_matches, best_matches;
        match_info      mi;
        
        buildCollisionAvoidingRandomTransform( mi.random_transform );
        transformModulePoles( mi.random_transform, transformed_M_poles );
        // convert this call to be inside StatefulEngine, in order to avoid passing the first two parameters
        match_module_to_host( *this->m, this->tree, transformed_M_poles, M_to_H );
        
        for( auto pole_and_vertex : M_to_H )
        {
            // cout << pole_and_vertex.first << ", " << pole_and_vertex.second << endl;
            current_matches.push_back( make_pair( pole_and_vertex.first, pole_and_vertex.second ));
        }
        EdgeCost c = get_best_subset( *this->m, current_matches, best_matches, no_glueings );
        graph_print( cout, c ) << endl;;
        mi.cost     = c;
        mi.matches  = best_matches;
        proposed_matches.push_back( mi );
        graph_print( cout, proposed_matches[i].cost ) << endl;;
    }
    
    std::cout << "after filling " << endl;
    assert( proposed_matches.size() > 0 );
    
    // find the best solution between the proposed ones
    size_t      selected = 0;
    EdgeCost    max_cost = proposed_matches[0].cost;
    graph_print( cout, max_cost ) << endl;
    for( int i = 1; i < proposed_matches.size(); ++i )
    {
        graph_print( cout, proposed_matches[i].cost ) << endl;;
        if( proposed_matches[i].cost < max_cost )
        {
            selected = i;
            max_cost = proposed_matches[i].cost;
        }
    }
    
    best_match.setMatchInfo( proposed_matches[selected] );
    cout << best_match.getMatchInfo().random_transform << endl;
    
    //    save_intermediate_result(host, TEST_PATH, 3);
    
    /* How this should work
     -) build a kdtree of the host
     -) for N times
     -)  compute a random transformation of the module (that brings it outside the host)
     it only applies the transformation to its poles
     -)  match each pole of the module to its closest point on the host
     -)  use the complete graph utilities in order to get the subset ( with given cardinality )
     that is the best fit between module and host
     -)  collect the configuration and its COST
     -) choose the configuration with smallest cost
     -) use svd rigid motion utility to align the module to the host given the matching between
     the choosen poles and vertices
     -) add a pole to the host in those choosen vertices that are not poles
     -) glue each choosen module's pole to its correspondent host's pole.
     -) smooth the added skeleton ( this is something that needs a little more work )
     */
    //    result = std::move( proposed_matches[selected].matches );
}


void StatefulEngine::glueModuleToHost(){
    alignUsingBestMatch();
    actualGlueing();
}








