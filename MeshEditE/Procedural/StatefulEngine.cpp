//#define TRACE
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
#include "MeshEditE/Procedural/Helpers/Plane.h"
#include "collision_detection.h"

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
    float   rand_max =  static_cast<float>( randomizer.max( )); //RAND_MAX / ( M_PI * 2.0 );
    float   x1       = static_cast<float>( randomizer( )) / rand_max,
            x2       = static_cast<float>( randomizer( )) / rand_max,
            x3       = static_cast<float>( randomizer( )) / rand_max;
    
    last_x1 = x1;   last_x2 = x2;   last_x3 = x3;
    
    t = ModuleAlignment::random_rotation_matrix_arvo( x1, x2, x3 );
}

void StatefulEngine::buildMainStructureKdTree(){
    assert( this->m != NULL );
    assert( this->tree == NULL );
    
    this->tree = new kD_Tree();
    ModuleAlignment::build_manifold_kdtree( (*this->m), mainStructure->getFreePoleSet(), *this->tree );
    
    for( auto& vid : mainStructure->getFreePoleSet( )){
        assert( is_pole( *m, vid ));
    }
    treeIsValid = true;
}


/********** MATCHING **********/

void StatefulEngine::matchModuleToHost( Module &candidate, VertexMatchMap& M_pole_to_H_vertex ){

    typedef vector<IdDistPair >                                Near_Pole_Vector;
    typedef map< VertexID, Near_Pole_Vector >                  Candidate_Neighbors;
    
    VertexSet           assigned_candidates, unassigned_poles;
    Candidate_Neighbors candidateNeighbors;
    // from module's pole to host's candidate
    VertexMatchMap      internal_match;
    
    // find the nearest host candidate for each module's pole
    // and, for each candidate matched, store its matched pole and distance
    for( const auto& id_and_info : candidate.getPoleInfoMap())
    {
        double      distance    = numeric_limits<double>::max();
        Vec3d       foundPos;
        VertexID    foundID = InvalidVertexID, second_choiceID = InvalidVertexID;
        bool        have_found  = (*tree).closest_point( id_and_info.second.geometry.pos , distance, foundPos, foundID );
        
        assert( have_found );
        assert( foundID != InvalidVertexID );
        assert( mainStructure->getFreePoleSet().count(foundID) > 0 );
        
        if( have_found ){
            const PoleInfo& host_pole_info  = mainStructure->getPoleInfo( foundID );
            const Vec3d&    n_candidate     = host_pole_info.geometry.normal;
            size_t          valence         = host_pole_info.geometry.valence;
            
            // add only if match is valid
#ifdef TRACE
            cout << "M : " << id_and_info.first << ")" << id_and_info.second.geometry.pos << "  #  " << id_and_info.second.geometry.normal << endl;
            cout << "H : " << foundID << ")" << foundPos << "  #  " << n_candidate << endl << endl;
            bool opposite_debug = opposite_directions( id_and_info.second.geometry.normal, n_candidate );
            cout << " opposite? " <<opposite_debug << "valence H : " << valence << " # M : " << id_and_info.second.geometry.valence<<endl;
#endif

            if( Module::poleCanMatch(host_pole_info, id_and_info.second)){
                if( opposite_directions( id_and_info.second.geometry.normal, n_candidate )){
                    
                    // instantiate vector if putting the first value
                    if( candidateNeighbors.count( foundID ) == 0 ){
                        candidateNeighbors[foundID] = Near_Pole_Vector();
                    }
                    
                    candidateNeighbors[foundID].push_back( make_pair( id_and_info.first, distance ));
                    assigned_candidates.insert( foundID );
                    // maps from module to main structure
                    internal_match[id_and_info.first] = foundID;
                }
            }
        }
    }
    
    // here I should process all the candidateNeighbors for which the number of matched
    // poles > 1 and try to assign them to not assigned candidates
    for( auto& item : candidateNeighbors ){
        assert( item.second.size() > 0 );
        if( item.second.size() == 1 ) {
            M_pole_to_H_vertex[item.second.front().first] = item.first;
        }
        else{
            // get the nearest matched pole
            IdDistPair& nearest = item.second.front();
            for( IdDistPair& id_dist : item.second ){
                if( IdDistPair_comparer( id_dist, nearest )){
                    nearest = id_dist;
                }
            }
            // save the nearest match into the output map
            M_pole_to_H_vertex[nearest.first] = item.first;
            assigned_candidates.insert( item.first );
//          put all the other poles that have found a neighbor into the unassigned set
            for( IdDistPair& id_dist : item.second ){
                if( id_dist.first != nearest.first ){ unassigned_poles.insert( id_dist.first ); }
            }
        }
    }
    
    // for each unassigned pole try to find a secondary nearest match
    for( VertexID unassigned : unassigned_poles ){
                
        assert( candidate.getPoleInfoMap().count( unassigned ) > 0 );
        const PoleGeometryInfo& pgi = candidate.getPoleInfo( unassigned ).geometry;

        VertexID second_cloesest = InvalidVertexID;
        
        if( findSecondClosest( unassigned, pgi, internal_match[unassigned], second_cloesest, assigned_candidates )){
            assert( mainStructure->getFreePoleSet().count(second_cloesest) > 0 );
            M_pole_to_H_vertex[unassigned] = second_cloesest;
        }
    }
    
    // sanity check
    map<VertexID, int> _candidates, _poles;
    for( const auto& item : M_pole_to_H_vertex ){
        _poles[item.first]          = 0;
        _candidates[item.second]    = 0;
    }
    
    for( auto& item : M_pole_to_H_vertex ){
        _poles[item.first] = _poles[item.first] + 1;
        _candidates[item.second] = _candidates[item.second] + 1;
        // each pole MUST be assigned to only one candidate and each candidate to only one pole
        assert( _poles.count( item.first ) <= 1);
        assert( _candidates.count( item.second ) <= 1);
//        assert( is_pole( *(candidateModule->m), item.first ));
        assert( candidateModule->isPole( item.first ));
        assert( is_pole( *m, item.second ));
#ifdef TRACE
        cout << item.first << " # " << item.second << endl;
#endif
    }
    
#ifdef TRACE
    cout << "candidates size " <<  _candidates.size() << endl;
#endif
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
#warning consider using squared length
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
//        Vec3d n_candidate = vertex_normal( *m, it->first );
//        n_candidate.normalize();
        Vec3d n_candidate = mainStructure->getPoleInfo(it->first).geometry.normal;
        
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



void StatefulEngine::buildOptimalAlignmentTransform( const Module& module, const std::vector<Match>& matches, CGLA::Mat4x4d& T ){
    Mat4x4d r, t;
    vector< Vec3d>     host_pos, module_pos, host_normals, module_normals;
    vector<double>     weights( matches.size(), 1.0 );
    for( auto& match : matches ){
        module_pos.push_back(       module.getPoleInfo( match.first ).geometry.pos );
        module_normals.push_back(   module.getPoleInfo(match.first).geometry.normal );
        host_pos.push_back(         mainStructure->getPoleInfo( match.second ).geometry.pos );
        host_normals.push_back(     mainStructure->getPoleInfo(match.second).geometry.normal );
        
        const Vec3d& mn = module.getPoleInfo(match.first).geometry.normal;
        const Vec3d& hn = mainStructure->getPoleInfo(match.second).geometry.normal;
    }
    
    if ( matches.size() > 1 ){
//        svd_rigid_motion( m   odule_pos, host_pos, r, t );
        svd_6d_rigid_motion( module_pos, host_pos, module_normals, host_normals, r, t );

        T = t * r;
    }
    else{
        Vec3d tr = host_pos.front() - module_pos.front();
        T = translation_Mat4x4d( tr );
    }

    truncateMat4x4d( T );
    checkMat4( T );
}

void StatefulEngine::buildNormalsAlignmentTransform( const Module& module, const std::vector<Match>& matches,
                                                    CGLA::Mat4x4d& T ){
    CGLA::Vec3d host_vec( 0 ), module_vec( 0 ), centroid( 0 );
    double      _;
    for( Match match : matches )
    {
        cout <<"centroid :" << centroid <<endl;
        centroid = centroid + module.getPoleInfo( match.first ).geometry.pos;
        Vec3d mn = module.getPoleInfo( match.first ).geometry.normal;
        Vec3d hn = mainStructure->getPoleInfo( match.second ).geometry.normal;
        
        assert( !( isnan( mn[0] ) || isnan( mn[1] ) || isnan( mn[2] )));
        assert( !( isnan( hn[0] ) || isnan( hn[1] ) || isnan( hn[2] )));
        
        checkVec3( centroid );
        
        module_vec  += mn;
        host_vec    += hn;
    }
    
    truncateVec3d( module_vec );
    truncateVec3d( host_vec );

    
#ifdef TRACE
    cout << "module vector :" << endl << module_vec << endl;
    cout << "host vector   :" << endl << module_vec << endl;
#endif
    
//    module_vec.normalize();
//    host_vec.normalize();
    centroid /= static_cast<double>( matches.size());
    
#ifdef TRACE
    cout << "module vector normalized :" << endl << module_vec << endl;
    cout << "host vector normalized   :" << endl << module_vec << endl;
    cout << "centroid                 :" << endl << centroid << endl;
#endif
    // need to carefully choose which centroid I should use.
    
    Mat4x4d tr_origin = translation_Mat4x4d( -centroid ),
            tr_back   = translation_Mat4x4d( centroid );
    Mat4x4d t_align   = alt_get_alignment_for_2_vectors( module_vec, host_vec );
            T         = tr_back * t_align * tr_origin;
    
    checkMat4( T );
    truncateMat4x4d( T );
    
#ifdef TRACE
    cout << "Best Match Normal Alignment " << endl << T << endl;
#endif

}

void StatefulEngine::applyRandomTransform(){
    //    cout << "transforming using : " << endl << best_match.getMatchInfo().random_transform << endl;
    Module &tm = candidateModule->getTransformedModule( best_match.getMatchInfo().random_transform);
    candidateModule = &tm;
    
    for( VertexID v : M_vertices ){
        m->pos( v) = best_match.getMatchInfo().random_transform.mul_3D_point( m->pos( v ));
    }
    
    assert( candidateModule->poleList.size() > 0 );
}


void StatefulEngine::applyOptimalAlignment(  ){
    Mat4x4d t;
    
    buildOptimalAlignmentTransform( *candidateModule, best_match.getMatchInfo().matches, t );
    
    
#ifdef TRACE
    cout << "Best Match Optimal (SVD) Alignment " << endl << t << endl;
#endif
    
    Module &tm = candidateModule->getTransformedModule( t );
    candidateModule = &tm;
    for( VertexID v : M_vertices ){
        m->pos( v) = t.mul_3D_point( m->pos( v ));
    }
    assert( candidateModule->poleList.size() > 0 );
}


void StatefulEngine::alignModuleNormalsToHost(){

    Mat4x4d t;
    buildNormalsAlignmentTransform( *candidateModule, best_match.getMatchInfo().matches, t );
    
    Module& tm = candidateModule->getTransformedModule( t );
    candidateModule = &tm;
    for( VertexID v : M_vertices ){
        m->pos( v) = t.mul_3D_point( m->pos( v ));
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
    vector<Match> remapped_matches;
    add_manifold(*m, *candidateModule->m, host_remap, module_remap, M_vertices );
    mainStructure->reAlignIDs( host_remap );
    candidateModule->reAlignIDs( module_remap );
    // remap matches ids
    for( Match& match : best_match.getMatchInfo().matches){

#ifdef TRACE
        cout << "old match : ( " << match.first << ", " << match.second << " )" << endl;
#endif

        assert( module_remap[match.first]  != InvalidVertexID );
        assert( host_remap[match.second]   != InvalidVertexID );

#ifdef TRACE
        cout << "new match : ( " << module_remap[match.first] << ", " << host_remap[match.second] << " )" << endl;
#endif
        remapped_matches.push_back( make_pair( module_remap[match.first], host_remap[match.second]) );
    }
    
    best_match.getMatchInfo().matches = std::move( remapped_matches );

//#ifdef TRACE
    for( Match& match : best_match.getMatchInfo().matches){
        cout << "match : ( " << candidateModule->getPoleInfo( match.first ).original_id << ", "
             << mainStructure->getPoleInfo( match.second ).original_id << " )" << endl;
    }
//#endif
    
    applyRandomTransform();
    candidateModule->updateDirections( *m );
    
    applyOptimalAlignment( );
    candidateModule->updateDirections( *m );

    
//    alignModuleNormalsToHost();
//    candidateModule->updateDirections( *m );
    
    // glue_matches
    Helpers::ModuleAlignment::glue_matches( *m, best_match.getMatchInfo().matches );
    
    // need to reset bounding sphere of the main structure
    Vec3d bsphere_center;
    double bsphere_radius;
    bsphere( *m, bsphere_center, bsphere_radius );
    mainStructure->setBoundingSphere( bsphere_center, bsphere_radius );
    mainStructure->glueModule( *m, *candidateModule, best_match.getMatchInfo().matches );
    mainStructure->updatePoleInfo( *m );
    
    IDRemap glue_remap;
    
    m->cleanup( glue_remap );
    mainStructure->reAlignIDs( glue_remap.vmap );

#ifdef TRACE
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
#endif
    consolidate();
    debugColorization();
}


void StatefulEngine::actualGlueing(){
    if( best_match.IsValid( )){
        
        // add manifold module to host
        set<VertexID> _h;
        IDRemap remap;
        VertexIDRemap host_remap, module_remap;
        vector<Match> remapped_matches;
        add_manifold(*m, *candidateModule->m, host_remap, module_remap, M_vertices );
        mainStructure->reAlignIDs( host_remap );
        candidateModule->reAlignIDs( module_remap );
        // remap matches ids
        for( Match& match : best_match.getMatchInfo().matches){
            remapped_matches.push_back( make_pair( module_remap[match.first], host_remap[match.second]) );
            
#ifdef TRACE
            cout << "old match : ( " << match.first << ", " << match.second << " )" << endl;
            cout << "new match : ( " << module_remap[match.first] << ", " << host_remap[match.second] << " )" << endl;
#endif
        }
        // glue_matches
        Helpers::ModuleAlignment::glue_matches( *m, remapped_matches );
        // mainStructure->glueModule
        mainStructure->glueModule( *m, *candidateModule, remapped_matches );
        
        IDRemap glue_remap;
        
        m->cleanup( glue_remap );
        mainStructure->reAlignIDs( glue_remap.vmap );
        
        consolidate();
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
    this->m = &host;
    this->mainStructure = new MainStructure();
    
    Module *starter = new Module( *m, 0 );
    std::vector<Procedural::Match> matches;
    mainStructure->glueModule( *m, *starter, matches);
    
    debugColorization();
}

void StatefulEngine::setHostFromModule( Procedural::Module &starter, HMesh::Manifold &host ){
    this->m = &host;
    this->mainStructure = new MainStructure();
    
    std::vector<Procedural::Match> matches;
    
    Vec3d bsphere_center;
    double bsphere_radius;
    bsphere( *m, bsphere_center, bsphere_radius );
    mainStructure->setBoundingSphere( bsphere_center, bsphere_radius );
    
    mainStructure->glueModule( *m, starter, matches);
    debugColorization();
}

void StatefulEngine::setModule(Procedural::Module &module){
    assert( module.m != NULL );
    this->candidateModule = &module;
    buildMainStructureKdTree();
}


void StatefulEngine::consolidate(){
    assert( this->m != NULL );
    assert( this->candidateModule->getPoleInfoMap().size() > 0 );
    // in this way you lose any reference to which vertices are from host and module

    candidateModule = NULL;
    transformedModules.clear();
    
    M_vertices.clear();
    treeIsValid = false;
    tree = NULL;
}



bool lessThan( double lhs, double rhs, int no_decimals, bool acceptEqual = false ){
    int _m = 10;
    // multiply for one time more in order to get an extra decimal and use that for comparison
    for( int i = 0; i <= no_decimals; ++i ){ _m *= 10; }

    double mult = static_cast<double>( _m );
    long eps = 5;

    double t_lhs = trunc( mult * lhs );
    double t_rhs = trunc( mult * rhs );
    long ti_lhs = static_cast<long>( t_lhs );
    long ti_rhs = static_cast<long>( t_rhs );
    long diff = ti_lhs - ti_rhs;
    // if diff < -eps is less
    // if diff > eps is greater
    // if -eps <= diff <= eps is equal
    
    cout <<
            lhs     << ", " << rhs      << endl <<
            t_lhs   << ", " << t_rhs    << endl <<
            ti_lhs  << ", " << ti_rhs   << endl <<
            diff    << endl;
    
    if( !acceptEqual ){
        return ( diff < -eps );
    }
    else{
        return diff < abs( eps );
    }
    
}

// this should maximize the number of matches while minimizing the total cost
// should add a little penalty to the single matches since their cost is 0
// can I minimize that in a least square sense?
size_t StatefulEngine::chooseBestFitting( const vector< match_info > proposed_matches, const vector< ExtendedCost > extendedCosts ) const{
    assert( proposed_matches.size() > 0 );
    size_t selected = 0;
    double best_sum_cost  = numeric_limits<double>::max();
    double d_penalty      = 0.0001;
    double current_cost   = 0.0;
    
    for( int i = 0; i < proposed_matches.size(); ++i )
    {
        ExtendedCost e = extendedCosts[i];
        size_t match_valency = proposed_matches[i].matches.size();
        if( match_valency == 1 ){
            current_cost = d_penalty;
        }else{
            double divider = static_cast<double>( match_valency );
            double squared_divider = divider * divider;            
            current_cost = ( e.first +  e.second.first+ e.second.second ) / squared_divider  ;
            
            cout << i << " # " << match_valency << " # " << divider << endl;
            
            cout << i << ") #" << match_valency << " => ( " << e.first
            << ", ( " << e.second.first << ", " << e.second.second << " ) = " << current_cost  << endl;

        }

        // clamp values and create INT for comparison
        if( lessThan( current_cost, best_sum_cost, 5, false ) ||
            ( lessThan( current_cost, best_sum_cost, 5, true ) &&
            ( proposed_matches[i].matches.size() > proposed_matches[selected].matches.size( )))){
        
            cout << "wins current "  << endl;
            best_sum_cost = current_cost;
            selected = i;
        }
        else{
            cout << "wins old best "  << endl;
        }
    }
    cout << "selected : " << selected << endl;
    return selected;
}


bool StatefulEngine::testMultipleTransformations(){
        assert( this->m != NULL );
    
    vector< matchesAndCost >                matches_vector;
    vector< Mat4x4d >                       Ts;
    vector< pair< size_t, ExtendedCost >>   stats;
    map<size_t, CandidateSubsetInfo>        best_subsets_for_valence;
    size_t                                  max_valence = 0;

    
    for( size_t d = 0; d < candidateModule->poleList.size(); ++d){ stats.push_back( make_pair( 0, make_pair( 0.0, make_pair( 0.0, 0.0 ))));}
    
    buildTransformationList( Ts );
//    return;
    
    assert( candidateModule->getPoleInfoMap().size() > 0 );
    
    for( int i = 0; i < Ts.size(); ++i ){

        VertexMatchMap  M_to_H;
        vector<Match>   current_matches, best_matches;
        
        assert( !isnan( Ts[i][1][1] )); // Need to discover why putting this assert prevents to assign NaN to mi.random_transform
        
        matchModuleToHost( transformedModules[i], M_to_H );
        
        if( M_to_H.size() <= 0 ) { continue; }
        
        stats[M_to_H.size() - 1].first++;
        
        double distance_sum = 0.0;
        
        for( auto& pole_and_vertex : M_to_H )
        {
            current_matches.push_back( make_pair( pole_and_vertex.first, pole_and_vertex.second ));
        }
        
        std::vector< SubsetResult > results;
        EdgeCost treshold = make_pair( 100.5, 100.5 );

        get_subsets( *mainStructure, transformedModules[i], current_matches, results, treshold );
        
//        if( results.size() == 0 ){ /* std::cout << "result set is empty" << endl; */ continue; }

        // debug!!!!
        if( results.size() < M_to_H.size() ){ /* std::cout << "result set is empty" << endl; */ continue; }

        // LEGACY
        assert( results.size() == 1 || results.front().matches.size() > results.back().matches.size( ));

        std:: cout << "|************  " << i << "************|" << endl;
        
        int count = 0;
        // BUILD SubsetInfo and check if they are "better matches" than the previous
        for( const SubsetResult& item : results ){
            CandidateSubsetInfo info;
            info.subset_match = item.matches;
            assert( info.subset_match.size() == item.matches.size( ));
            info.T_random           = Ts[i];
            
            //            Module& step1 = transformedModules[i];
            //            buildOptimalAlignmentTransform( step1, info.subset_match, info.T_align );
            //
            //            Module& step2 = step1.getTransformedModule( info.T_align );
            //            buildNormalsAlignmentTransform( step2, info.subset_match, info.T_normals );
            
            Module& step1 = transformedModules[i];
            buildOptimalAlignmentTransform( step1, info.subset_match, info.T_align );
            
            mat4Copy( info.T_align * info.T_random, info.T_complete );

//            mat4Copy( info.T_normals * info.T_align * info.T_random, info.T_complete );

            checkMat4( info.T_complete );
            truncateMat4x4d( info.T_complete );

//            info.transformed_module = step2.getTransformedModule( info.T_normals );

            info.transformed_module = step1.getTransformedModule( info.T_align );

            // calculate costs
            cout <<  "\t *** subset " << count++ << endl;
            for( Match m : item.matches ){
                
                const PoleInfo& module_pole_info = info.transformed_module.getPoleInfo( m.first );
                const PoleInfo& host_pole_info   = mainStructure->getPoleInfo( m.second );
                
                info.module_ball_radius = candidateModule->bsphere_radius;
                info.distance_cost +=
                    ( module_pole_info.geometry.pos - host_pole_info.geometry.pos ).length();
                
                info.distance_normal_angle +=
                    normal_distance( dot( module_pole_info.geometry.normal, host_pole_info.geometry.normal ));
                
                // calculation of the cost relative to the alignmento of the preferred direction
                // needs to be explained
                if( module_pole_info.anisotropy.is_defined && host_pole_info.anisotropy.is_defined ){
                    double dot_value = dot( module_pole_info.anisotropy.direction, host_pole_info.anisotropy.direction );
                    cout << module_pole_info.anisotropy.direction << endl << host_pole_info.anisotropy.direction << endl << endl;
                    
                    // if anisotropy is bilateral for both I need to consider as Optimal also if directions are opposite
                    if( module_pole_info.anisotropy.is_bilateral && host_pole_info.anisotropy.is_bilateral ){
                        info.distance_alignment += bilateral_anisotropy_distance( dot_value );
                    }else{
                        // if it works
                        info.distance_alignment += anisotropy_distance( dot_value );
                    }
                }
            }
            
//            if( info.distance_alignment > 0.0005 ){ continue; }
            
            info.calculateTotalCost( true );

            subsetsInfo.push_back( info );
            
            cout << count << ") " << endl
                 << info.distance_cost << ", "
                << info.distance_alignment << ", "
                << info.distance_normal_angle << endl;
            
            
            // save the less expensive subsets clustered by valence.
            if( best_subsets_for_valence.count( info.matchValence() ) == 0 ){
                best_subsets_for_valence[info.matchValence()] = info;
            }
            else{
//                    if( info.getTotalCost() < best_subsets_for_valence[info.matchValence()].getTotalCost()){
                if( CandidateSubsetInfo::betterThen( info, best_subsets_for_valence[info.matchValence()] )){
                    best_subsets_for_valence[info.matchValence()] = info;                        
                }
            }
            
            if( info.matchValence() > max_valence ) { max_valence = info.matchValence(); }
            
            cout
            << "match with valence : " << info.matchValence()
            << " with cost : ( " << info.distance_cost
            << ", " << info.distance_normal_angle
            << ", " << info.distance_alignment << " ) => "
            << info.getTotalCost( )  << endl << endl
            << "****************************************" << endl   ;
            
        }
    }
//    assert( best_subsets_for_valence.count( max_valence ) > 0 );
    if( best_subsets_for_valence.size() == 0  ){
        consolidate();
        return false;
    }
    CandidateSubsetInfo best_subset = best_subsets_for_valence[max_valence];
    
    for( int i = max_valence; i > 0; --i ){
        if( best_subsets_for_valence.count( i ) > 0 ){
            
            cout << "match with valence : " << best_subsets_for_valence[i].matchValence()
            << " with cost : ( " << best_subsets_for_valence[i].distance_cost
            << ", " << best_subsets_for_valence[i].distance_normal_angle
            << ", " << best_subsets_for_valence[i].distance_alignment << " ) => "
            << best_subsets_for_valence[i].getTotalCost()      << endl
            << "while best was : ( "   << best_subset.matchValence()
            << ", "                    << best_subset.getTotalCost()
            << " )" << endl;
            
            if( CandidateSubsetInfo::betterThen( best_subsets_for_valence[i], best_subset )){
//            if( best_subsets_for_valence[i].getTotalCost( ) < best_subset.getTotalCost( )){
                best_subset = best_subsets_for_valence[i];
            }
        }
    }
    match_info      mi;
//    mi.random_transform = best_subset.T_random;
    mat4Copy( best_subset.T_random, mi.random_transform );
    mi.matches = std::move( best_subset.subset_match );
    best_match.setMatchInfo( mi );
    
    
    if( !best_match.IsValid() ){
        cout << "unable to find a feasible solution of " << transformedModules.size() << " transformed modules " << endl;
        consolidate();
        return false;
    }
    

#ifdef TRACE
    cout << "Best Match random transform: " << selected << endl
         << best_match.getMatchInfo().random_transform
         << " Proposed Matches " << endl;
    for( Match& match :  best_match.getMatchInfo().matches ){
        cout << match.first << " # " << match.second << endl;
        Vec3d normal1 = vertex_normal( *m, match.first );
        Vec3d normal2 = vertex_normal( *m, match.second );
        normal1.normalize();
        normal2.normalize();
        cout << "1) " << normal1 << endl << "2) " << normal2 << endl;
        double angle = get_angle(normal1, normal2);
        cout << " angle : " << angle << " and opposite? :   " << opposite_directions( normal1, normal2 ) << endl;

    }
#endif
    
    // debug
    for( size_t d = 0; d < candidateModule->poleList.size(); ++d){
        std::cout << d+1 << ") " << stats[d].first << endl;
    }

    return true;
}


void StatefulEngine::glueModuleToHost(){
    alignUsingBestMatch();
    actualGlueing();
}

void StatefulEngine::buildTransformationList( vector< Mat4x4d> &transformations ){

    size_t skipped = 0;
    
#ifdef TRACE
    cout << "Building transformations set " << endl;
#endif
    
    size_t no_M_poles = candidateModule->getPoleInfoMap().size();
    transformations.clear();
    
    const vector<VertexID> &_candidates = mainStructure->getFreePoles();
    
    size_t no_candidates = _candidates.size();
    size_t no_m_poles    = candidateModule->poleList.size();
    size_t H_starter     = randomizer() % no_candidates;
    Mat4x4d t_origin     = translation_Mat4x4d( - candidateModule->bsphere_center );
    
    for( int i = 0; i < no_candidates; ++i ){
        
        size_t actual_i = ( i + H_starter ) % no_candidates;
        VertexID H_pole = _candidates[ actual_i ];
//        VertexID H_pole = _candidates[ 0 ];
        
        const PoleInfo& H_pole_info = mainStructure->getPoleInfo( H_pole );
        
        size_t M_starter = randomizer() % no_M_poles;
        // MODULE POLES LOOP
        for( int j = 0; j < no_M_poles; ++j ){
            
            Mat4x4d T;
            size_t actual_j = ( j + M_starter ) % no_m_poles;
            VertexID M_pole = candidateModule->poleList[ actual_j ];
//            VertexID M_pole = candidateModule->poleList[ 0 ];
            
            assert( candidateModule->getPoleInfoMap().count( M_pole ) > 0 );
            const PoleInfo& pinfo  = candidateModule->getPoleInfo( M_pole );
            

            if( !( Module::poleCanMatch( pinfo, H_pole_info ))){
                skipped += 16;
                continue;
            }
            
            if( pinfo.original_id == H_pole_info.original_id ){
                cout << "check this";
            }
            
            cout << " testing  : ( " << pinfo.original_id << ", " << H_pole_info.original_id << " )" << endl;
            
            // this can be changed to use the new 6d svd rotation
            // align normals
            Mat4x4d t_align = alt_get_alignment_for_2_vectors( pinfo.geometry.normal, H_pole_info.geometry.normal );
            
            Vec3d m_pole_step1 = ( t_align * t_origin ).mul_3D_point( pinfo.geometry.pos );
            
            // generate K rotations along the candidate normal, according to the case
            // standard case, at least one pole is unconstrained
            int no_steps = 16;
            double step = M_PI_4 / 2.0;
            double curr_angle = ( randomizer() % 16 ) * step;
//            double curr_angle = 0;
            
            Vec3d m_pole_anis_dir_aligned;
            // both of them
            if( pinfo.anisotropy.is_defined && H_pole_info.anisotropy.is_defined ){

                Vec3d moved_pole     = ( t_align * t_origin ).mul_3D_point( pinfo.geometry.pos );
                Vec3d to_H_pole      = H_pole_info.geometry.pos - moved_pole;
                Mat4x4d tr_to_H_pole = translation_Mat4x4d( to_H_pole );
                
                Mat4x4d Tdir = ( tr_to_H_pole * t_align * t_origin );
                truncateMat4x4d( Tdir );
                
                        m_pole_anis_dir_aligned = mul_3D_dir( Tdir, pinfo.anisotropy.direction );
                        moved_pole              = Tdir.mul_3D_point( pinfo.geometry.pos );
                Vec3d   moved_normal            = mul_3D_dir( Tdir, pinfo.geometry.normal );
                
                truncateVec3d( m_pole_anis_dir_aligned );
                
                Plane alpha( H_pole_info.geometry.pos, H_pole_info.geometry.normal );
                cout << "before projecting on plane " << m_pole_anis_dir_aligned << endl;
                m_pole_anis_dir_aligned = alpha.projectDirection( m_pole_anis_dir_aligned );
                cout << "after projecting on plane " << m_pole_anis_dir_aligned << endl;
                
                double  anis_angle              = get_angle( m_pole_anis_dir_aligned, H_pole_info.anisotropy.direction );
                
                
                cout << "applyng " << endl
                << Tdir << endl
                << " to "<< endl
                << pinfo.anisotropy.direction << endl
                << "results in " << m_pole_anis_dir_aligned << endl;

                //test truncate anis_angle
                anis_angle = truncateDouble3( anis_angle );
                
                if( pinfo.anisotropy.is_bilateral && H_pole_info.anisotropy.is_bilateral ){
                    no_steps    = 2;
                    step        = M_PI;
                    curr_angle = ( randomizer() % 2 ) == 0 ? anis_angle : anis_angle + step;
                }
                else{ // only one is bilateral
                    no_steps    = 1;
                    step        = 2.0 * M_PI;
                    curr_angle  = anis_angle;
                }
                

                // debug
                assert(  dot( pinfo.anisotropy.direction, pinfo.geometry.normal ) - 1.0 < 0.0000001 );
                assert(  dot( H_pole_info.anisotropy.direction, H_pole_info.geometry.normal ) - 1.0 < 0.0000001 );
//                Plane alpha( H_pole_info.geometry.pos, H_pole_info.geometry.normal );
                Plane beta ( moved_pole, moved_normal );
                cout << " plane alpha - main structure " << endl << alpha.toString() << endl;
                cout << " plane beta  - main structure " << endl << beta.toString() << endl;
//                assert( alpha.OnPlane( H_pole_info.geometry.pos ));
//                assert( alpha.OnPlane( moved_pole ));
//                assert( beta.OnPlane( moved_pole + m_pole_anis_dir_aligned ));
//                assert( alpha.OnPlane( H_pole_info.geometry.pos + H_pole_info.anisotropy.direction ));
                
                cout
                << "current angle " << curr_angle << endl
                << "before rotation around host axis " << endl
                << pinfo.anisotropy.direction << pinfo.anisotropy.direction.length() << endl
                << m_pole_anis_dir_aligned << m_pole_anis_dir_aligned.length() << endl
                << H_pole_info.anisotropy.direction << H_pole_info.anisotropy.direction.length() << endl;
            }
//            curr_angle = -curr_angle;
            

            // build and save rotations
            for ( int i = 0; i < no_steps; ++i, curr_angle +=step ) {
                Mat4x4d rot = get_rotation_mat4d( H_pole_info.geometry.normal, curr_angle );
                truncateMat4x4d( rot );
                Vec3d rot_result = mul_3D_dir( rot, m_pole_anis_dir_aligned );
                double rot_result_dot = dot( rot_result, H_pole_info.anisotropy.direction );
                if( anisotropy_distance( rot_result_dot ) > ARITH_EPS ){
                    rot = get_rotation_mat4d( H_pole_info.geometry.normal, -curr_angle );
                }

                Vec3d m_pole_step2 = rot.mul_3D_point( m_pole_step1 );
                
                Vec3d   to_H_pole = H_pole_info.geometry.pos - m_pole_step2;
                Mat4x4d tr_to_H_pole = translation_Mat4x4d( to_H_pole );
#ifdef TRACE
                cout << "rotation with axis " << H_pole_normal << " and angle : " << curr_angle << rot << endl;
                cout << "Translate to pole " << endl << tr_to_H_pole << endl;
                cout << "pole translated aligned rotated : " << m_pole_step2 << endl;
                cout << "to_h_pole " << to_H_pole << endl
                << "translation mat4 " << endl << tr_to_H_pole;
#endif
                
                truncateMat4x4d( tr_to_H_pole );
                truncateMat4x4d( t_align );
                truncateMat4x4d( tr_to_H_pole );
                truncateMat4x4d( t_origin );
                
                T = tr_to_H_pole * rot * t_align * t_origin;

                truncateMat4x4d( T );

                checkMat4( t_origin );
                checkMat4( t_align );
                checkMat4( rot );
                checkMat4( tr_to_H_pole );
                checkMat4( T );

#ifdef TRACE
                cout << "complete transform" << endl <<  T;
#endif
                
                Module& t_module = this->candidateModule->getTransformedModule( T );
                // end - debug
                double anis_dot_after_rotation = dot( t_module.getPoleInfo(M_pole).anisotropy.direction, H_pole_info.anisotropy.direction );
                if( fabs( 1.0 - anis_dot_after_rotation ) > 0.03 ){
                    cout << "pssss ehi check here! angle is : "
                         << anis_dot_after_rotation << " # " << acos( anis_dot_after_rotation ) << endl
                         << " direction history " << endl;
                    Vec3d step1 = mul_3D_dir( t_origin, pinfo.anisotropy.direction );
                    Vec3d step2 = mul_3D_dir( t_align, step1 );
                    Vec3d step3 = mul_3D_dir( rot, step2);
                    Vec3d step4 = mul_3D_dir( tr_to_H_pole, step3);
                    
                    cout
                        << pinfo.anisotropy.direction << endl
                        << step1 << endl
                        << step2 << endl
                        << step3 << endl
                        << step4 << endl
                        << anis_dot_after_rotation;
                }
                
                // skip if there is a collision.
                // need to improve it
                if( mainStructure->isColliding( t_module )){
                    cout << " skipping configuration - COLLISION " << endl;
                    continue;
                }
                
                cout << "after rotation around host axis " << endl
                << t_module.getPoleInfo(M_pole).anisotropy.direction
                << t_module.getPoleInfo(M_pole).anisotropy.direction.length() << endl
                << H_pole_info.anisotropy.direction << H_pole_info.anisotropy.direction.length() << endl;


                transformations.push_back( T );
                transformedModules.push_back( t_module );
            }
        }
    }
    
    assert( transformations.size() == transformedModules.size( ));
    cout << endl << transformations.size() << " configurations generated and " << skipped << " skipped" << endl;
}

size_t StatefulEngine::noFreePoles(){
    return this->mainStructure->getFreePoles().size();
}

void StatefulEngine::debugColorization(){
    
    Vec3f v_color( 0.0, 0.0, 1.0 );
    Vec3f f_color( 0.9, 0.9, 0.9 );
    Vec3f e_color( 0.0, 0.0, 0.0 );
    Vec3f d_color( 1.0, 0.0, 0.0 ); // directions
    Vec3f p_color( 0.0, 1.0, 1.0 ); // pole
    for( auto vit = m->vertices().begin(); vit != m->vertices().end(); ++vit )
    {
        assert( m->in_use( *vit ));
        bool is_pole = mainStructure->getFreePoleSet().count( *vit ) > 0;
        
        
        GLGraphics::DebugRenderer::vertex_colors[*vit] = is_pole ? p_color : v_color;
        
        Walker w = m->walker(*vit);
        for(; !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            if( is_pole && mainStructure->getPoleInfo(*vit).anisotropy.is_defined ){
                if( w.vertex() == mainStructure->getPoleInfo(*vit).anisotropy.directionID ){
                    GLGraphics::DebugRenderer::edge_colors[w.halfedge()]        = d_color;
                    GLGraphics::DebugRenderer::edge_colors[w.opp().halfedge()]  = d_color;
                }
                else{
                    GLGraphics::DebugRenderer::edge_colors[w.halfedge()]        = e_color;
                    GLGraphics::DebugRenderer::edge_colors[w.opp().halfedge()]  = e_color;
                }
            }
        }
    }
    for( auto fit = m->faces().begin(); fit != m->faces().end(); ++fit ){
        GLGraphics::DebugRenderer::face_colors[*fit] = f_color;
    }
}