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
    this->m         = NULL;
    this->tree      = NULL;
    this->module    = NULL;
    unsigned seed   = chrono::system_clock::now().time_since_epoch().count();
    randomizer.seed( seed );
    treeIsValid     = false;
    dim_constraint  = DimensionalityConstraint::Constrained_3D;
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
    double  dir_len = dir.length();
    double  length  = H_radius + M_radius - dir.length();
    
    cout << "Dir : " << dir << " # with length : " << dir_len << endl;
    
    if( fabs( dir_len ) < 0.0000000001 ){
        dir = Vec3d( last_x1, last_x2, last_x3 );
    }
    
    dir.normalize();
    dir *= length;
    
    tr = CGLA::translation_Mat4x4d( dir );
}


void StatefulEngine::buildCollisionAvoidingRandomTransform( CGLA::Mat4x4d &t ){
    Mat4x4d rot, tr;
    // build random transform
    buildRandomRotation( rot );
    cout << "Random Rotation" << endl << rot ;
    // build collision avoiding translation
    buildCollisionAvoidingTranslation( rot, tr );
    cout << "Random Translation" << endl << tr ;
    // return composition
    t = tr * rot;
}


void StatefulEngine::buildHostKdTree(){
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    assert( this->tree == NULL );
    
//    ModuleAlignment::build_manifold_kdtree( (*this->m), this->H_vertices, this->tree );
    this->tree = new kD_Tree();
    ModuleAlignment::build_manifold_kdtree( (*this->m), this->H_candidates.getCandidates(), *this->tree );
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


/********** MATCHING **********/

void StatefulEngine::matchModuleToHost( Procedural::PoleInfoMap& poleInfoMap, VertexMatchMap& M_pole_to_H_vertex ){

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
        assert( H_candidates.getCandidates().count(foundID) > 0 );
        
        if( have_found ){
            Vec3d n_candidate = vertex_normal( *m, foundID );
            n_candidate.normalize();
            size_t valence = Geometry::valence( *m, foundID );
            
            // add only if match is valid
            
            cout << "M : " << id_and_info.first << ")" << id_and_info.second.geometry.pos << "  #  " << id_and_info.second.geometry.normal << endl;
            cout << "H : " << foundID << ")" << foundPos << "  #  " << n_candidate << endl << endl;
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
            assert( H_candidates.getCandidates().count(second_cloesest) > 0 );
            M_pole_to_H_vertex[unassigned] = second_cloesest;
        }
    }
    
    // sanity check
    map<VertexID, int> _candidates, _poles;
    for( auto item : M_pole_to_H_vertex ){
        _poles[item.first] = 0;
        _candidates[item.second] = 0;
    }
    
    for( auto item : M_pole_to_H_vertex ){
        _poles[item.first] = _poles[item.first] + 1;
        _candidates[item.second] = _candidates[item.second] + 1;
        // each pole MUST be assigned to only one candidate and each candidate to only one pole
        assert( _poles.count( item.first ) <= 1);
        assert( _candidates.count( item.second ) <= 1);
        assert( is_pole( *m, item.first ));
        assert( is_pole( *m, item.second ));
        cout << item.first << " # " << item.second << endl;
    }
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
    cout << "transforming using : " << endl << best_match.getMatchInfo().random_transform << endl;
    
    for( auto mv : M_vertices ) {
        m->pos( mv ) = best_match.getMatchInfo().random_transform.mul_3D_point( m->pos( mv ));
    }
}


void StatefulEngine::applyOptimalAlignment(){
    Mat4x4d R, T;
    vector< VertexID > host_v, module_p;
    for( auto m : best_match.getMatchInfo().matches )
    {
        module_p.push_back(m.first);
        host_v.push_back(m.second);
        Vec3f red( 1, 0, 0 ), blue(0, 1, 0);
        
        // use debug colors to exploit the matched vertices
        GLGraphics::DebugRenderer::vertex_colors[m.first] = red;
        GLGraphics::DebugRenderer::vertex_colors[m.second] = blue;
        
    }
    //#warning there is an active return here!
    //    return;
    
    //    save_intermediate_result(m, TEST_PATH , 1);
    svd_rigid_motion( *m, module_p, *m, host_v, R, T );
    Mat4x4d t = T * R;
    cout << "Best Match Optimal (SVD) Alignment " << endl << t << endl;
    for( auto mv : M_vertices ) { m->pos( mv ) = t.mul_3D_point( m->pos( mv )); }
}


void StatefulEngine::alignModuleNormalsToHost(){
    CGLA::Vec3d host_vec( 0 ), module_vec( 0 ), centroid( 0 );
    double      _;
    for( Match match : best_match.getMatchInfo().matches )
    {
        centroid += m->pos( match.first );
        Vec3d mn = vertex_normal( *m, match.first );
        Vec3d hn = vertex_normal( *m, match.second );
        mn.normalize();
        hn.normalize();
        module_vec  += mn;
        host_vec    -= hn;
    }
    
    cout << "module vector :" << endl << module_vec << endl;
    cout << "host vector   :" << endl << module_vec << endl;
    
    module_vec.normalize();
    host_vec.normalize();
    centroid /= best_match.getMatchInfo().matches.size();
    
    cout << "module vector normalized :" << endl << module_vec << endl;
    cout << "host vector normalized   :" << endl << module_vec << endl;
    cout << "centroid                 :" << endl << centroid << endl;

    // need to carefully choose which centroid I should use.
    //    bsphere( m, module_IDs, centroid, _ );
    Mat4x4d tr_origin = translation_Mat4x4d( -centroid ),
            tr_back   = translation_Mat4x4d( centroid );
    Mat4x4d t_align = get_alignment_for_2_vectors( module_vec, host_vec );
    Mat4x4d t = tr_back * t_align * tr_origin;
    cout << "Best Match Normal Alignment " << endl << t << endl;
    for( VertexID v : M_vertices )
    {
        m->pos(v) = t.mul_3D_point( m->pos( v ));
    }
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


void StatefulEngine::actualGlueing(){
    if( best_match.IsValid( )){
        Helpers::ModuleAlignment::glue_matches( *m, best_match.getMatchInfo().matches );
        consolidate();
    }
}


void StatefulEngine::fillCandidateSet(){
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    
    for( VertexID vid : H_vertices ){
        if( is_pole( *m, vid )) {
            CandidateInfo _;
            H_candidates.insert( vid, _ );
        }
    }
}


/// If some of the selected vertices from the host are not poles, then they must be converted into poles.
/// Since doing this changes the structure of the mesh, we need to clean up the Manifold
/// and update the sets of vertices' IDs of the host and the module
void StatefulEngine::addNecessaryPoles( )
{
    vector<VertexID>                    host_poles_to_remap, host_poles_remapped;
    set<VertexID>                       module_vs_remapped,  host_vs_remapped;
    HMesh::VertexAttributeVector<int>   vertex_selection;
    IDRemap                             remap;
    vector<Match>& matches = best_match.getMatchInfo().matches;
    // aggiungo i poli e salvo il loro id
    for( int i = 0; i < matches.size(); ++i )
    {
        // get the valence of the poles
        
        if( !is_pole( *m, matches[i].second ))
        {
            // clean the selection
            for( auto v : m->vertices( )) { vertex_selection[v] = 0; }
            // calculate the right size accordingly to the valence of the corresponding pole

            int size = current_glueing_target / 4 ;
            assert(( current_glueing_target % 4 ) == 0 );
            //            int size = -1;
            //            assert( size != -1 );
            VertexID poleID = add_branch( *m, matches[i].second, size, vertex_selection );
            assert( poleID != InvalidVertexID );
            matches[i].second = poleID;

        }
    }
    m->cleanup( remap );
    // re-align the saved IDs
    for( auto p : remap.vmap )
    {
        if( p.second == InvalidVertexID ) { continue; }
        if( M_vertices.count( p.first ) > 0 )    { module_vs_remapped.insert( p.second ); }
        else                                     { host_vs_remapped.insert( p.second );   }
    }
    
    // change
    for( int i = 0; i < matches.size(); ++i )
    {
        matches[i].first  = remap.vmap[matches[i].first];
        matches[i].second = remap.vmap[matches[i].second];
    }
    //    selected    = std::move( host_poles_remapped );
    M_vertices   = std::move( module_vs_remapped );
    H_vertices   = std::move( host_vs_remapped );
    
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
    edge_info.Update( m );
}

void StatefulEngine::setModule( Manifold &module ){
    assert( this->m != NULL );
    assert( this->H_vertices.size() == 0 );
    assert( this->M_vertices.size() == 0 );
    add_manifold( (*this->m), module, H_vertices, M_vertices );
    fillCandidateSet();
    buildHostKdTree();

    // LOAD MODULE INFO
    this->module = new Module( "", 0 );
    Vec3d centroid;
    double radius;
    bsphere( *m, M_vertices, centroid, radius );
    this->module->bsphere_center = centroid;
    this->module->bsphere_radius = radius;
    
    for( VertexID vid : M_vertices ){
        if( is_pole( *m, vid )){
            PoleInfo pi;
            pi.geometry.valence = valence( *m, vid );
            pi.geometry.pos     = m->pos( vid );
            Vec3d n = vertex_normal( *m, vid );
            n.normalize();
            pi.geometry.normal = n;
            this->module->poleList.push_back( vid );
            this->module->poleInfoMap[vid] = pi;
        }
    }
}

void StatefulEngine::consolidate(){
    assert( this->m != NULL );
    assert( this->H_vertices.size() > 0 );
    assert( this->M_vertices.size() > 0 );
    assert( this->module->poleInfoMap.size() > 0 );
    // in this way you lose any reference to which vertices are from host and module
    module = NULL;
    transformedModules.clear();
    
    H_vertices.clear();
    M_vertices.clear();
    H_candidates.clear();
    treeIsValid = false;
    edge_info.Invalidate();

    kD_Tree* temp = tree;
    tree = NULL;
    delete temp;
}


void StatefulEngine::testMultipleTransformations( int no_tests, size_t no_glueings ){
        assert( this->m != NULL );
        assert( this->H_vertices.size() > 0 );
        assert( this->M_vertices.size() > 0 );
        assert( no_tests > 0 );
    
    current_glueing_target = no_glueings;
    
    vector< matchesAndCost >    matches_vector;
    vector< match_info >        proposed_matches;
    vector< ExtendedCost >      ExtendetCosts;
    vector< Mat4x4d >           Ts;
    
    buildTransformationList( Ts );
    return;
    
    assert( module->poleInfoMap.size() > 0 );
    
    for( int i = 0; i < Ts.size(); ++i ){

        VertexMatchMap  M_to_H;
        vector<Match>   current_matches, best_matches;
        match_info      mi;
        
        assert( !isnan( Ts[i][1][1] )); // Need to discover why putting this assert prevents to assign NaN to mi.random_transform
        
        mi.random_transform = Ts[i];
        assert( !isnan( mi.random_transform[1][1] ));

//        cout << " Ts[i] " << Ts[i] << endl;
//        cout << " Using as Random Transform" << endl
//             << mi.random_transform
//             << " test assignment : " << endl
//             << aggg;
//        
//        // can be removed when solved
//        cout << "poles with normals " << endl;
//        for( auto item : transformedModules[i].poleInfoMap )
//        {
//            cout << item.first << ")" << item.second.geometry.pos << "  #  " << item.second.geometry.normal << endl;
//        }
        
        matchModuleToHost( transformedModules[i].poleInfoMap, M_to_H );
        
        // if the number of matches is lower than the target, just skip this configuration
        if( M_to_H.size() < no_glueings ) { continue; }
        
        double distance_sum = 0.0;
        
        for( auto pole_and_vertex : M_to_H )
        {
            // cout << pole_and_vertex.first << ", " << pole_and_vertex.second << endl;
            current_matches.push_back( make_pair( pole_and_vertex.first, pole_and_vertex.second ));
        }

        EdgeCost c = get_best_subset( *this->m, current_matches, best_matches, no_glueings );
        
        for( auto match : best_matches ){
            distance_sum += ( transformedModules[i].poleInfoMap[match.first].geometry.pos - m->pos(match.second)).length();
        }

        mi.cost     = c;
        ExtendetCosts.push_back( make_pair( distance_sum, c ));

        cout << "configuration " << i << " has cost : "
        << mi.cost.first << ", " << mi.cost.second << ", " << distance_sum << endl;

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

    cout << "Building transformations set " << endl;
    
    size_t no_M_poles = module->poleInfoMap.size();
    Vec3d H_centroid, M_centroid;
    double H_radius, M_radius;
    transformations.clear();
    bsphere( *m, H_vertices, H_centroid, H_radius );
    bsphere( *m, M_vertices, M_centroid, M_radius );
    
//    cout << "module centroid :" << M_centroid << endl;
    
    Mat4x4d tr_to_origin    = translation_Mat4x4d( -M_centroid ),
            tr_to_centroid  = translation_Mat4x4d( M_centroid );
    
//    cout << "Translate to origin " << endl << tr_to_origin << endl << " and back " << endl << tr_to_centroid;
    
    vector<VertexID> _candidates = H_candidates.getCandidateVector();
    
    assert(_candidates.size() == H_candidates.getCandidates().size());
    
    size_t no_candidates = H_candidates.getCandidates().size();
    size_t no_m_poles    = module->poleList.size();
    size_t H_starter = randomizer() % no_candidates;
    
//    for( VertexID H_pole : H_candidates.getCandidates() ){
    for( int i = 0; i < _candidates.size(); ++i ){
        
        size_t actual_i = ( i + H_starter ) % no_candidates;
//        VertexID H_pole = _candidates[ actual_i ];
        VertexID H_pole = _candidates[ 5 ];

        Mat4x4d T;
        Vec3d   H_pole_pos      = m->pos( H_pole );
        Vec3d   H_pole_normal   = vertex_normal( *m, H_pole );
        H_pole_normal.normalize();
        
        size_t M_starter = randomizer() % no_M_poles;

        for( int j = 0; j < _candidates.size(); ++j ){
            
            size_t actual_j = ( j + M_starter ) % no_m_poles;
//            VertexID M_pole = module->poleList[ actual_j ];
              VertexID M_pole = module->poleList[ 0 ];
            
            cout << " polo : " << M_pole << endl;

            assert( module->poleInfoMap.count(M_pole) > 0 );
            PoleInfo pinfo  = module->poleInfoMap[M_pole];
            
            // align normals
//            Mat4x4d t_align = get_alignment_for_2_vectors( pinfo.geometry.normal, H_pole_normal );
            Mat4x4d t_align = alt_get_alignment_for_2_vectors( pinfo.geometry.normal, H_pole_normal );
            cout << "M_pole # pos : " << pinfo.geometry.pos << " # normal : " << pinfo.geometry.normal << endl
                 << "H_pole # pos : " << H_pole_pos         << " # normal : " << H_pole_normal << endl;
            
            cout << "alignment" << t_align << endl;
            
            
            // bring the transformed M_pole to H_pole
            Vec3d   t_align_M_pole = t_align.mul_3D_point( pinfo.geometry.pos );
            Vec3d   to_H_pole = H_pole_pos - t_align_M_pole;
//            Vec3d   to_H_pole = H_pole_pos - pinfo.geometry.pos;
            
            Mat4x4d tr_to_H_pole = translation_Mat4x4d( to_H_pole );
//            cout << "Translate to pole " << endl << tr_to_H_pole << endl;
            
            cout << "t_aligned M_pole # pos : " << t_align_M_pole << endl;
            cout << "to_h_pole " << to_H_pole << endl
                 << "translation mat4 " << endl << tr_to_H_pole;

            
            // generate K rotations along the candidate normal
            double step = M_PI_4 / 2.0;
            double curr_angle = 0;
            // build and save rotations
            for ( int i = 0; i < 1; ++i, curr_angle +=step ) {
                Mat4x4d rot = get_rotation_mat4d( H_pole_normal, curr_angle );
                
//                cout << "rotation with axis " << H_pole_normal << " and angle : " << curr_angle << rot << endl;
                
                Vec3d centroid = t_align.mul_3D_point( M_centroid );
                Mat4x4d tr_to_origin    = translation_Mat4x4d( -H_pole_pos ),
                        tr_to_centroid  = translation_Mat4x4d( H_pole_pos );


//                T = tr_to_H_pole * tr_to_origin * rot * tr_to_centroid * t_align;
//                T=  tr_to_H_pole * tr_to_centroid * t_align * tr_to_origin;
//                T = t_align;
                T = tr_to_H_pole * t_align;
                cout << "complete transform" << endl <<  T;

//                cout << "complete transformation" << T << endl;
                
                transformations.push_back( T );
                assert( !isnan( T[1][1] ));
                
                Module t_module( "", 0 );
                
                // SAVE THE MODULE TRANSFORMATION
                for( auto item : module->poleInfoMap ){
                    PoleGeometryInfo pgi = item.second.geometry;
                    PoleInfo t_pi;
                    t_pi.geometry.valence = pgi.valence;
                    t_pi.geometry.pos     = T.mul_3D_point( pgi.pos );
                    
                    // TODO vector must be multiplied expanding to 4d vector with w=0!!!!!!
                    Vec4d normal( pgi.normal[0], pgi.normal[1], pgi.normal[2], 0.0 );
                    Vec4d t_normal4 = T * normal;
                    Vec3d t_normal( t_normal4[0], t_normal4[1], t_normal4[2] );
                    t_normal.normalize();
//                  cout << pgi.normal << " becomes : " << t_normal << endl;
                    t_pi.geometry.normal  = t_normal;
                    t_module.poleInfoMap[ item.first ] = t_pi;
                }
                assert( t_module.poleInfoMap.size() == module->poleInfoMap.size( ));
                transformedModules.push_back( t_module );
            }
        }
    }
    
    for( auto mv : M_vertices ) {
        m->pos( mv ) = transformations[8].mul_3D_point( m->pos( mv ));
        if( is_pole(*m, mv)){ cout << mv << ") " << m->pos(mv) << endl; }
    }
    
    assert( transformations.size() == transformedModules.size( ));
}