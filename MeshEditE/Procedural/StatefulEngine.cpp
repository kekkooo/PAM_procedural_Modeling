    //
//  StatefulEngine.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 17/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "StatefulEngine.h"

#include <GEL/GLGraphics/ManifoldRenderer.h>

#include "polarize.h"

#include "MeshEditE/Procedural/Helpers/manifold_copy.h"
#include "MeshEditE/Procedural/Helpers/geometric_properties.h"
#include "MeshEditE/Procedural/Helpers/svd_alignment.h"
#include "MeshEditE/Procedural/Operations/structural_operations.h"

#include "Test.h"

using namespace std;

using namespace HMesh;
using namespace CGLA;

using namespace Procedural::Engines;
using namespace Procedural::Helpers;
using namespace Procedural::Geometry;
using namespace Procedural::Helpers::ModuleAlignment;
using namespace Procedural::GraphMatch;
using namespace Procedural::Operations::Structural;
using namespace Procedural::EngineHelpers;

/*=========================================================================*
 *                     PRIVATE FUNCTIONS                                   *
 *=========================================================================*/

StatefulEngine::StatefulEngine()
{
    this->m                 = NULL;
    this->tree              = NULL;
    this->candidateModule   = NULL;
    this->mainStructure     = NULL;
    unsigned seed   = chrono::system_clock::now().time_since_epoch().count();
    randomizer.seed( seed );
    treeIsValid     = false;

}

/********** UTILITIS **********/

/// returns true if L < R
bool IdDistPair_comparer( IdDistPair &l, IdDistPair &r ) { return ( l.second < r.second ); }


void StatefulEngine::buildRandomRotation( CGLA::Mat4x4d &t ){
    // must be initialized
    float rand_max =  static_cast<float>( randomizer.max( )); //RAND_MAX / ( M_PI * 2.0 );
    float x1 = static_cast<float>( randomizer( )) / rand_max,
    x2 = static_cast<float>( randomizer( )) / rand_max,
    x3 = static_cast<float>( randomizer( )) / rand_max;
    
    last_x1 = x1;   last_x2 = x2;   last_x3 = x3;
    
    t = ModuleAlignment::random_rotation_matrix_arvo( x1, x2, x3 );
}



void StatefulEngine::buildMainStructureKdTree(){
    assert( this->m != NULL );
    assert( this->tree == NULL );
    
//    ModuleAlignment::build_manifold_kdtree( (*this->m), this->H_vertices, this->tree );
    this->tree = new kD_Tree();
    ModuleAlignment::build_manifold_kdtree( (*this->m), mainStructure->getFreePoleSet(), *this->tree );
    for( auto vid : mainStructure->getFreePoleSet( )){
        assert( is_pole(*m, vid));
    }
    treeIsValid = true;
}


/********** MATCHING **********/

void StatefulEngine::matchModuleToHost( const Procedural::PoleInfoMap& poleInfoMap, VertexMatchMap& M_pole_to_H_vertex ){

    typedef vector<IdDistPair >                                Near_Pole_Vector;
    typedef map< VertexID, Near_Pole_Vector >                  Candidate_Neighbors;
    
    VertexSet           assigned_candidates, unassigned_poles;
    Candidate_Neighbors candidateNeighbors;
    // from module's pole to host's candidate
    VertexMatchMap      internal_match;
    
    // find the nearest host candidate for each module's pole
    // and, for each candidate matched, store its matched pole and distance
    for( auto id_and_info : poleInfoMap)
    {
        double      distance    = numeric_limits<double>::max();
        Vec3d       foundPos;
        VertexID    foundID = InvalidVertexID, second_choiceID = InvalidVertexID;
        bool        have_found  = (*tree).closest_point( id_and_info.second.geometry.pos , distance, foundPos, foundID );
        
        assert( have_found );
        assert( foundID != InvalidVertexID );
        assert( mainStructure->getFreePoleSet().count(foundID) > 0 );
        
        if( have_found ){
            cout << " HAVE FOUND "<< endl;
            Vec3d n_candidate = vertex_normal( *m, foundID );
            n_candidate.normalize();
            size_t valence = Geometry::valence( *m, foundID );
            
            // add only if match is valid
            
            cout << "M : " << id_and_info.first << ")" << id_and_info.second.geometry.pos << "  #  " << id_and_info.second.geometry.normal << endl;
            cout << "H : " << foundID << ")" << foundPos << "  #  " << n_candidate << endl << endl;
            bool opposite_debug = opposite_directions( id_and_info.second.geometry.normal, n_candidate );
            cout << " opposite? " <<opposite_debug << "valence H : " << valence << " # M : " << id_and_info.second.geometry.valence<<endl;
            if( opposite_directions( id_and_info.second.geometry.normal, n_candidate )
             && valence  == id_and_info.second.geometry.valence ){
                
                // instantiate vector if putting the first value
                if( candidateNeighbors.count( foundID ) == 0 ){
                    candidateNeighbors[foundID] = Near_Pole_Vector();
                }
                candidateNeighbors[foundID].push_back( make_pair( id_and_info.first, distance ));
                assigned_candidates.insert( foundID );
                internal_match[id_and_info.first] = foundID;
            }
        }
    }
    
    // here I should process all the candidateNeighbors for which the number of matched
    // poles > 1 and try to assign them to not assigned candidates
    for( auto item : candidateNeighbors ){
        assert( item.second.size() > 0 );
        if( item.second.size() == 1 ) {
            M_pole_to_H_vertex[item.second.front().first] = item.first;
        }
        else{
            // get the nearest matched pole
            IdDistPair nearest = item.second.front();
            for( IdDistPair id_dist : item.second ){
                if( IdDistPair_comparer( id_dist, nearest )){
                    nearest = id_dist;
                }
            }
            // save the nearest match into the output map
            M_pole_to_H_vertex[nearest.first] = item.first;
            assigned_candidates.insert( item.first );
            // put all the other poles into the unassigned set
            for( IdDistPair id_dist : item.second ){
                if( id_dist.first != nearest.first ){ unassigned_poles.insert( id_dist.first ); }
            }
        }
    }
    
    // for each unassigned pole try to find a secondary nearest match
    for( VertexID unassigned : unassigned_poles ){
        assert(poleInfoMap.count(unassigned) > 0 );
        PoleGeometryInfo pgi = poleInfoMap[unassigned].geometry;

        VertexID second_cloesest = InvalidVertexID;
        if( findSecondClosest( unassigned, pgi, internal_match[unassigned], second_cloesest, assigned_candidates )){
            assert( mainStructure->getFreePoleSet().count(second_cloesest) > 0 );
            M_pole_to_H_vertex[unassigned] = second_cloesest;
        }
    }
    
    // sanity check
    map<VertexID, int> _candidates, _poles;
    for( auto item : M_pole_to_H_vertex ){
        _poles[item.first]          = 0;
        _candidates[item.second]    = 0;
    }
    
    for( auto item : M_pole_to_H_vertex ){
        _poles[item.first] = _poles[item.first] + 1;
        _candidates[item.second] = _candidates[item.second] + 1;
        // each pole MUST be assigned to only one candidate and each candidate to only one pole
        assert( _poles.count( item.first ) <= 1);
        assert( _candidates.count( item.second ) <= 1);
        assert( is_pole( *(candidateModule->m), item.first ));
        assert( is_pole( *m, item.second ));
        cout << item.first << " # " << item.second << endl;
    }
    
    cout << "candidates size " <<  _candidates.size() << endl;
}


bool StatefulEngine::findSecondClosest( const VertexID &pole, const PoleGeometryInfo &pgi, const VertexID &closest, VertexID &second_closest, VertexSet &assigned ){
    assert( assigned.count( closest ) > 0 );
    // consider closest
    Vec3d       closest_pos = m->pos( closest );
    // find the mean distance from 1-neighbors
    double      mean_dist   = numeric_limits<double>::max();
    size_t      valence     = 0;
    for( Walker w = m->walker( closest ); !w.full_circle(); w = w.circulate_vertex_ccw())
    {
        double dist = ( closest_pos - m->pos( w.vertex( ))).length( );
        mean_dist += dist;
        ++valence;
    }
    // take that distance as radius and use it with tree.in_sphere with center in closest
    double              query_radius = mean_dist / (float)valence;
    vector<VertexID>    in_sphere_ID;
    vector<Vec3d>       in_sphere_points;
    IDsDistsVector      ids_and_dists;
    
    (*tree).in_sphere( closest_pos, query_radius, in_sphere_points, in_sphere_ID );
    // take the nearest point
    assert( in_sphere_ID.size() == in_sphere_points.size());
    
    for( int i = 0; i < in_sphere_points.size(); ++i){
        ids_and_dists.push_back( make_pair( in_sphere_ID[i], (closest_pos - in_sphere_points[i] ).length( )));
    }
    
    // sort in ascending order accordingly to the distance
    std::sort( ids_and_dists.begin(), ids_and_dists.end(), IdDistPair_comparer );
    IDsDistsVector::iterator it      = ids_and_dists.begin();
    bool                    done    = false;
    while ( !done && it != ids_and_dists.end( ) )
    {        
        // check if normals are compatible
        Vec3d n_candidate = vertex_normal( *m, it->first );
        n_candidate.normalize();

        done    = (( assigned.count( it->first ) == 0 )
                  && opposite_directions( pgi.normal, n_candidate ))
                  && Geometry::valence( *m, it->first ) == pgi.valence;
        if( !done ) { ++it; }
    }
    if( done ){
        second_closest = it->first;
        assigned.insert( it->first );
    }
    return done;
}


/********** APPLICATION OF TRANSFORMATIONS **********/

void StatefulEngine::applyRandomTransform(){
//    cout << "transforming using : " << endl << best_match.getMatchInfo().random_transform << endl;
    Module &tm = candidateModule->getTransformedModule( best_match.getMatchInfo().random_transform);
    candidateModule = &tm;
    
    for( VertexID v : M_vertices ){
        m->pos( v) = best_match.getMatchInfo().random_transform.mul_3D_point( m->pos( v ));
    }
    
    assert( candidateModule->poleList.size() > 0 );
}


void StatefulEngine::applyOptimalAlignment(){
    Mat4x4d R, T;
    vector< VertexID > host_v, module_p;
    for( auto m : best_match.getMatchInfo().matches ){
        module_p.push_back(m.first);
        host_v.push_back(m.second);
    }
    //#warning there is an active return here!
    //    return;
    
    //    save_intermediate_result(m, TEST_PATH , 1);
    svd_rigid_motion( *m, module_p, *(candidateModule->m), host_v, R, T );
    Mat4x4d t = T * R;
//    cout << "Best Match Optimal (SVD) Alignment " << endl << t << endl;
    
    Module &tm = candidateModule->getTransformedModule( t );
    candidateModule = &tm;
    for( VertexID v : M_vertices ){
        m->pos( v) = best_match.getMatchInfo().random_transform.mul_3D_point( m->pos( v ));
    }
    assert( candidateModule->poleList.size() > 0 );
}


void StatefulEngine::alignModuleNormalsToHost(){
    CGLA::Vec3d host_vec( 0 ), module_vec( 0 ), centroid( 0 );
    double      _;
    for( Match match : best_match.getMatchInfo().matches )
    {
        centroid += candidateModule->m->pos( match.first );
        Vec3d mn = candidateModule->getPoleInfo(match.first).geometry.normal;
        Vec3d hn = vertex_normal( *m, match.second );
        mn.normalize();
        hn.normalize();
        module_vec  += mn;
        host_vec    += hn;
    }
    
//    cout << "module vector :" << endl << module_vec << endl;
//    cout << "host vector   :" << endl << module_vec << endl;
//    
    module_vec.normalize();
    host_vec.normalize();
    centroid /= best_match.getMatchInfo().matches.size();
    
//    cout << "module vector normalized :" << endl << module_vec << endl;
//    cout << "host vector normalized   :" << endl << module_vec << endl;
//    cout << "centroid                 :" << endl << centroid << endl;

    // need to carefully choose which centroid I should use.
    //    bsphere( m, module_IDs, centroid, _ );
    Mat4x4d tr_origin = translation_Mat4x4d( -centroid ),
            tr_back   = translation_Mat4x4d( centroid );
    Mat4x4d t_align = alt_get_alignment_for_2_vectors( module_vec, host_vec );
    Mat4x4d t = tr_back * t_align * tr_origin;
//    cout << "Best Match Normal Alignment " << endl << t << endl;
    
    Module& tm = candidateModule->getTransformedModule( t );
    candidateModule = &tm;
    for( VertexID v : M_vertices ){
        m->pos( v) = best_match.getMatchInfo().random_transform.mul_3D_point( m->pos( v ));
    }

    assert( candidateModule->poleList.size() > 0 );
}


void StatefulEngine::alignUsingBestMatch( ){
    if( best_match.IsValid( )){
        applyRandomTransform();
        applyOptimalAlignment( );
        alignModuleNormalsToHost();
    }
    else{
        cout << "Warning : Best Match is not valid " << endl;
    }
}



void StatefulEngine::glueCurrent(){
    // add manifold module to host
    set<VertexID> _h;
    IDRemap remap;
    VertexIDRemap host_remap, module_remap;
    vector<GraphMatch::Match> remapped_matches;
    add_manifold(*m, *candidateModule->m, host_remap, module_remap, M_vertices );
    mainStructure->reAlignIDs( host_remap );
    candidateModule->reAlignIDs( module_remap );
    // remap matches ids
    for( Match match : best_match.getMatchInfo().matches){
        remapped_matches.push_back( make_pair( module_remap[match.first], host_remap[match.second]) );
    }
    
    best_match.getMatchInfo().matches = std::move( remapped_matches );
    
    applyRandomTransform();
    applyOptimalAlignment();
    alignModuleNormalsToHost();
    
    // glue_matches
    Helpers::ModuleAlignment::glue_matches( *m, best_match.getMatchInfo().matches );
    // mainStructure->glueModule
    mainStructure->glueModule( *candidateModule, best_match.getMatchInfo().matches );
    
    IDRemap glue_remap;
    
    m->cleanup( glue_remap );
    mainStructure->reAlignIDs( glue_remap.vmap );
    
    for( VertexID v : (*candidateModule).poleList ){
        if( mainStructure->getFreePoleSet().count(v) > 0 ){
            cout << "on module" << endl;
            cout << candidateModule->getPoleInfo(v).geometry.pos << endl;
            cout << candidateModule->getPoleInfo(v).geometry.normal << endl;
            cout << "on main structure" << endl;
            cout << m->pos(v) << endl;
            Vec3d normal = vertex_normal(*m, v);
            normal.normalize();
            cout << normal << endl;
            cout << endl<<endl;
            
        }
    }
    consolidate();
}






void StatefulEngine::actualGlueing(){
    if( best_match.IsValid( )){
        
        // add manifold module to host
        set<VertexID> _h;
        IDRemap remap;
        VertexIDRemap host_remap, module_remap;
        vector<GraphMatch::Match> remapped_matches;
        add_manifold(*m, *candidateModule->m, host_remap, module_remap, M_vertices );
        mainStructure->reAlignIDs( host_remap );
        candidateModule->reAlignIDs( module_remap );
        // remap matches ids
        for( Match match : best_match.getMatchInfo().matches){
            remapped_matches.push_back( make_pair( module_remap[match.first], host_remap[match.second]) );
            cout << "old match : ( " << match.first << ", " << match.second << " )" << endl;
            cout << "new match : ( " << module_remap[match.first] << ", " << host_remap[match.second] << " )" << endl;
        }
        // glue_matches
        Helpers::ModuleAlignment::glue_matches( *m, remapped_matches );
        // mainStructure->glueModule
        mainStructure->glueModule( *candidateModule, remapped_matches );
        
        IDRemap glue_remap;
        
        m->cleanup( glue_remap );
        mainStructure->reAlignIDs( glue_remap.vmap );
        
        for( VertexID v : (*candidateModule).poleList ){
            if( mainStructure->getFreePoleSet().count(v) > 0 ){
                cout << "on module" << endl;
                cout << candidateModule->getPoleInfo(v).geometry.pos << endl;
                cout << candidateModule->getPoleInfo(v).geometry.normal << endl;
                cout << "on main structure" << endl;
                cout << m->pos(v) << endl;
                Vec3d normal = vertex_normal(*m, v);
                normal.normalize();
                cout << normal << endl;
                cout << endl<<endl;

            }
        }
        consolidate();
    }
}

/// If some of the selected vertices from the host are not poles, then they must be converted into poles.

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
    this->mainStructure = new MainStructure();
    
    Module *starter = new Module( *m, 0 );
    std::vector<Procedural::GraphMatch::Match> matches;
    mainStructure->glueModule( *starter, matches);
    
    edge_info.Update( m );
}


void StatefulEngine::setModule(Procedural::Module &module){
    assert( module.m != NULL );
    this->candidateModule = &module;
    buildMainStructureKdTree();
}

void StatefulEngine::setModule( Manifold &module ){
//    assert( this->m != NULL );
//    assert( this->M_vertices.size() == 0 );
//    IDRemap remapper;
//    add_manifold( (*this->m), module, H_vertices, M_vertices, remapper );
//    mainStructure->reAlignIDs( remapper );
//    buildMainStructureKdTree();
//
//    // LOAD MODULE INFO
//    this->candidateModule = new Module();
//    Vec3d centroid;
//    double radius;
//    bsphere( *m, M_vertices, centroid, radius );
//    this->candidateModule->bsphere_center = centroid;
//    this->candidateModule->bsphere_radius = radius;
//    
//    for( VertexID vid : M_vertices ){
//        if( is_pole( *m, vid )){
//            PoleInfo pi;
//            pi.geometry.valence = valence( *m, vid );
//            pi.geometry.pos     = m->pos( vid );
//            Vec3d n             = vertex_normal( *m, vid );
//            n.normalize();
//            pi.geometry.normal  = n;
//            this->candidateModule->poleList.push_back( vid );
//            this->candidateModule->poleInfoMap[vid] = pi;
//        }
//    }
}

void StatefulEngine::consolidate(){
    assert( this->m != NULL );
    assert( this->M_vertices.size() > 0 );
    assert( this->candidateModule->getPoleInfoMap().size() > 0 );
    // in this way you lose any reference to which vertices are from host and module

    Module* tmp_module = candidateModule;
    candidateModule = NULL;
    transformedModules.clear();
    
    M_vertices.clear();
    treeIsValid = false;
    edge_info.Invalidate();

    kD_Tree* temp = tree;
    tree = NULL;
    delete temp;
    delete tmp_module;
}


void StatefulEngine::testMultipleTransformations( int no_tests, size_t no_glueings ){
        assert( this->m != NULL );
        assert( no_tests > 0 );
    
    current_glueing_target = no_glueings;
    
    vector< matchesAndCost >    matches_vector;
    vector< match_info >        proposed_matches;
    vector< ExtendedCost >      ExtendetCosts;
    vector< Mat4x4d >           Ts;
    
    buildTransformationList( Ts );
//    return;
    
    assert( candidateModule->getPoleInfoMap().size() > 0 );
    
    for( int i = 0; i < Ts.size(); ++i ){

        VertexMatchMap  M_to_H;
        vector<Match>   current_matches, best_matches;
        match_info      mi;
        
        assert( !isnan( Ts[i][1][1] )); // Need to discover why putting this assert prevents to assign NaN to mi.random_transform
        
        mi.random_transform = Ts[i];
        assert( !isnan( mi.random_transform[1][1] ));

        cout << " Ts[i] " << Ts[i] << endl;

        // can be removed when solved
        cout << "poles with normals " << endl;
        for( auto item : transformedModules[i].getPoleInfoMap() )
        {
            cout << item.first << ")" << item.second.geometry.pos << "  #  " << item.second.geometry.normal << endl;
        }
        
        matchModuleToHost( transformedModules[i].getPoleInfoMap(), M_to_H );
        
        // if the number of matches is lower than the target, just skip this configuration
        if( M_to_H.size() < no_glueings ) { continue; }
        
        double distance_sum = 0.0;
        
        for( auto pole_and_vertex : M_to_H )
        {
            cout << pole_and_vertex.first << ", " << pole_and_vertex.second << endl;
            current_matches.push_back( make_pair( pole_and_vertex.first, pole_and_vertex.second ));
        }

        EdgeCost c = get_best_subset( *this->m, current_matches, best_matches, no_glueings );
        
        for( auto match : best_matches ){
            distance_sum += ( transformedModules[i].getPoleInfo(match.first).geometry.pos - m->pos(match.second)).length();
        }

        mi.cost     = c;
        ExtendetCosts.push_back( make_pair( distance_sum, c ));

//        cout << "configuration " << i << " has cost : "
//        << mi.cost.first << ", " << mi.cost.second << ", " << distance_sum << endl;

        mi.matches  = std::move( best_matches );
        proposed_matches.push_back( std::move( mi ));
}
    

#warning this assert MUST be uncommented and this return deleted
//    return;
    assert( proposed_matches.size() > 0 );
    
    // find the best solution between the proposed ones
    size_t      selected = 0;
    ExtendedCost min_cost = ExtendetCosts[0];
    for( int i = 1; i < proposed_matches.size(); ++i )
    {
        if( ExtendetCosts[i] < min_cost )
        {
            selected = i;
            min_cost = ExtendetCosts[i];
        }
    }
    
    best_match.setMatchInfo( proposed_matches[selected] );

//    cout << "Best Match random transform: " << selected << endl
//         << best_match.getMatchInfo().random_transform
//         << " Proposed Matches " << endl;
//    for( Match match :  best_match.getMatchInfo().matches ){
//        cout << match.first << " # " << match.second << endl;
//        Vec3d normal1 = vertex_normal( *m, match.first );
//        Vec3d normal2 = vertex_normal( *m, match.second );
//        normal1.normalize();
//        normal2.normalize();
//        cout << "1) " << normal1 << endl << "2) " << normal2 << endl;
//        double angle = get_angle(normal1, normal2);
//        cout << " angle : " << angle << " and opposite? :   " << opposite_directions( normal1, normal2 ) << endl;
//
//}
    
    }


void StatefulEngine::glueModuleToHost(){
    alignUsingBestMatch();
    actualGlueing();
}

void StatefulEngine::buildTransformationList( vector< Mat4x4d> &transformations ){

//    cout << "Building transformations set " << endl;
    
    size_t no_M_poles = candidateModule->getPoleInfoMap().size();
    transformations.clear();
    
//    cout << "module centroid :" << M_centroid << endl;
    
//    cout << "Translate to origin " << endl << tr_to_origin << endl << " and back " << endl << tr_to_centroid;
    
    const vector<VertexID> &_candidates = mainStructure->getFreePoles();
    
    size_t no_candidates = _candidates.size();
    size_t no_m_poles    = candidateModule->poleList.size();
    size_t H_starter = randomizer() % no_candidates;
    Mat4x4d t_origin = translation_Mat4x4d( - candidateModule->bsphere_center );
    
    for( int i = 0; i < no_candidates; ++i ){
        
        size_t actual_i = ( i + H_starter ) % no_candidates;
        VertexID H_pole = _candidates[ actual_i ];
//        VertexID H_pole = _candidates[ 0 ];

        Vec3d   H_pole_pos      = m->pos( H_pole );
        Vec3d   H_pole_normal   = vertex_normal( *m, H_pole );
        H_pole_normal.normalize();
        
        size_t M_starter = randomizer() % no_M_poles;
        // MODULE POLES LOOP
        for( int j = 0; j < no_M_poles; ++j ){
            
            Mat4x4d T;
            size_t actual_j = ( j + M_starter ) % no_m_poles;
            VertexID M_pole = candidateModule->poleList[ actual_j ];
//            VertexID M_pole = candidateModule->poleList[ 0 ];
            
//            cout << " polo : " << M_pole << endl;

            assert( candidateModule->getPoleInfoMap().count(M_pole) > 0 );
            PoleInfo pinfo  = candidateModule->getPoleInfo(M_pole);
            
            // align normals
            Mat4x4d t_align = alt_get_alignment_for_2_vectors( pinfo.geometry.normal, H_pole_normal );
            
            Vec3d m_pole_step1 = ( t_align * t_origin ).mul_3D_point( pinfo.geometry.pos );

#warning randomize this. also the step
            // generate K rotations along the candidate normal
            double step = M_PI_4 / 2.0;
            double curr_angle = 0;
            // build and save rotations
            for ( int i = 0; i < 16; ++i, curr_angle +=step ) {
                Mat4x4d rot = get_rotation_mat4d( H_pole_normal, curr_angle );
                Vec3d m_pole_step2 = rot.mul_3D_point( m_pole_step1 );
                
                Vec3d   to_H_pole = H_pole_pos - m_pole_step2;
                Mat4x4d tr_to_H_pole = translation_Mat4x4d( to_H_pole );

//                cout << "rotation with axis " << H_pole_normal << " and angle : " << curr_angle << rot << endl;
//                cout << "Translate to pole " << endl << tr_to_H_pole << endl;
//                cout << "pole translated aligned rotated : " << m_pole_step2 << endl;
//                cout << "to_h_pole " << to_H_pole << endl
//                << "translation mat4 " << endl << tr_to_H_pole;

                T = tr_to_H_pole * rot * t_align * t_origin;
//                cout << "complete transform" << endl <<  T;
                
                transformations.push_back( T );
                assert( !isnan( T[1][1] ));
                
                Module t_module = this->candidateModule->getTransformedModule( T );
                
                transformedModules.push_back( t_module );
            }
        }
    }
    
//    for( auto mv : M_vertices ) {
//        m->pos( mv ) = transformations[8].mul_3D_point( m->pos( mv ));
//        if( is_pole(*m, mv)){ cout << mv << ") " << m->pos(mv) << endl; }
//    }
    
    assert( transformations.size() == transformedModules.size( ));
}