//
//  manifold_copy.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 06/11/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "manifold_copy.h"
#include <GEL/CGLA/Vec3d.h>

#include "polarize.h"

using namespace HMesh;
using namespace std;
using namespace CGLA;

namespace Procedural{
    namespace Helpers{

        
void add_manifold( Manifold &m, Manifold& other, VertexIDRemap &m_poles_remap,
                   VertexIDRemap &other_poles_remap, set<VertexID> & other_IDs_in_result){
    
    // save the set of the original vertices IDs of destination
    set<VertexID>                       m_IDs_set;
    set<VertexID>                       m_poles;
    map< FaceID, vector<FaceID> >       face_steps_map;
    map<FaceID, FaceID>                 old_to_new_faces;
    map<VertexID, vector<VertexID> >    possibleNewIDs;

    IDRemap remap;
    m_poles_remap.clear();
    other_poles_remap.clear();
    
    // save each vertexID of m into m_IDs_set, and its poles into m_poles
    for( VertexID v : m.vertices( )) {
        m_IDs_set.insert( v );
        if( is_pole( m, v )) {
            m_poles.insert( v );
            
            // save the pole onering along with it.
            Walker w = m.walker( v );
            while( !w.full_circle( )){
                m_poles.insert(w.vertex());
                w = w.circulate_vertex_ccw();
            }
        }
    }
    // for each face of other save the connnectivity and add it to m
    for( auto f : other.faces( ))
    {
        vector< Vec3d > vs;
        Walker          w       = other.walker(f);
        int             steps   = -1;
        VertexID        pole    = InvalidVertexID;
        
        // iterate along the face borders in order to know which are the neighbouring faces
        // and save them along with the number of steps needed to get there
        while( !w.full_circle( ))
        {
            if( face_steps_map.count( f ) == 0 ) { face_steps_map[f] = vector< FaceID >(); }
            face_steps_map[f].push_back( w.opp().face() );
            
            if( is_pole( other, w.prev().vertex( ))){
                steps = vs.size();
                pole  = w.prev().vertex();
            }
            
            // prendo il prev, perché altrimenti parto una rotazione avanti
            vs.push_back( other.pos( w.prev().vertex( )));
            w = w.next();
            
            // if vertex is pole, keep its step number
        }

        steps = vs.size() - steps - 1;
        
        // add the face keeping track of the number of m's vertices before and after the face addition
        size_t  pre_size    = m.no_vertices();
        FaceID  new_f       = m.add_face( vs );
        size_t  post_size   = m.no_vertices();
        
        assert( pre_size + vs.size() == post_size );
        
        if( pole != InvalidVertexID ){
            
            auto end = m.vertices_end();
            --end;

            for( int i = steps; i > 0; --i, --end);

            if( possibleNewIDs.count( pole ) == 0){
                possibleNewIDs[pole] = vector<VertexID>();
            }
            possibleNewIDs[pole].push_back( *end );
        }
        
        // map save pole ID into m after adding the face. use step number to know where to find it on m.vertices

        assert(new_f != InvalidFaceID );
        old_to_new_faces[f] = new_f;
    }
    // for each face of other - merge the faces using as reference the number of steps
    // needed to the walker to get the right edge
    for( auto f : other.faces( ))
    {
        FaceID  new_f       = old_to_new_faces[f];
        auto    neighbors   = face_steps_map[f];
        Walker  fw          = m.walker( new_f );
        
        for( int steps = 0; steps < neighbors.size(); ++steps )
        {
            // find the neighor face
            FaceID old_f_neighbor = neighbors[steps];
            FaceID new_neighbor   = old_to_new_faces[old_f_neighbor];
            
            // find N that is the index at face_steps_map[old_f_neighbor that gives f
            int i = 0;
            for(  ; i < neighbors.size() && face_steps_map[old_f_neighbor][i] != f; ++i ); // loops without doing anything
            assert( face_steps_map[old_f_neighbor][i] == f );
            
            Walker neighbor_w = m.walker(new_neighbor);
            for( int j = 0; j < i; ++j ) { neighbor_w = neighbor_w.next(); }
            m.stitch_boundary_edges( neighbor_w.opp().halfedge(), fw.opp().halfedge( ));
            fw = fw.next();
        }
    }
    
    // map from Other_ID to m_ID
    VertexIDRemap other_to_m_pre_cleanup;
    
    for( const auto& item : possibleNewIDs ){
        int sum = 0;
        for( VertexID candidate : item.second ){
            if ( m.in_use( candidate )){
                ++sum;
                other_to_m_pre_cleanup[item.first] = candidate;
            }
        }
        assert( sum == 1);
    }
    
    // test
    map<VertexID, VertexID> poleNeighborRemap;
    
    // for each pole of other copied to m
    for( auto item : other_to_m_pre_cleanup  ){
        Walker o_walker = other.walker( item.first );
        Walker m_walker = m.walker( item.second );
        
        Vec3d   m_starter_pos       = m.pos( m_walker.vertex() );
        bool    starter_not_found   = true;
        
        // loop until you find a correspondance between the neighbors of the pole on the two manifolds
        while( starter_not_found ){
            
            Vec3d   o_pos               = other.pos( o_walker.vertex( ));
                    starter_not_found   = ( o_pos - m_starter_pos ).length() > 0.00001;
            if ( starter_not_found ){
                o_walker = o_walker.circulate_vertex_ccw();
            }
        }
    
        Vec3d o_pos = other.pos( o_walker.vertex( ));
        
        assert(( o_pos - m_starter_pos ).length() < 0.00002 );
        assert( !o_walker.full_circle( ));
        
        // once the corresponding vertices have been found, iterate the pole one-ring
        // and save the remapping.
        
        while( !m_walker.full_circle()){
            
            poleNeighborRemap[o_walker.vertex()] = m_walker.vertex();
            
            m_walker = m_walker.circulate_vertex_ccw();
            o_walker = o_walker.circulate_vertex_ccw();
        }
    }
    // save it on the shared list.
    for( auto item : poleNeighborRemap  ){
        other_to_m_pre_cleanup[item.first] = item.second;
    }
    
    // endtest
    
    
    m.cleanup( remap );
    
    // find which occurence of new pole ID is the good one and put it into the other_remap
    
    // this remaps only the original vertices of m
    for( VertexID p : m_poles ){
        if( m_poles.count( p )) { m_poles_remap[p] = remap.vmap[p]; };
    }
    
    // remap both the origina vertices of m and its onering
    for( const auto& item : other_to_m_pre_cleanup ){
        VertexID newID = remap.vmap[ item.second ];
        assert( newID != InvalidVertexID );
        other_poles_remap[item.first] = newID;
    }

    other_IDs_in_result.clear();
    for( auto p : remap.vmap ){
        if( !m_IDs_set.count( p.first )) { other_IDs_in_result.insert( p.second ); }
    }

}
    


void add_manifold( Manifold &m, Manifold &other, set<VertexID> &m_IDs_in_result,
                  set<VertexID> & other_IDs_in_result, IDRemap &remap ){
    // save the set of the original vertices IDs of destination
    set<VertexID>                   m_IDs_set;
    map< FaceID, vector<FaceID> >   face_steps_map;
    map<FaceID, FaceID>             old_to_new_faces;
    m_IDs_in_result.clear();
    other_IDs_in_result.clear();
    
    for( VertexID v : m.vertices( )) { m_IDs_set.insert(v); }
    // for each face of source save the connnectivity and add it to destination
    for( auto f : other.faces( ))
    {
        vector< Vec3d > vs;
        Walker          w   = other.walker(f);
        while( !w.full_circle( ))
        {
            if( face_steps_map.count(f) == 0 ) { face_steps_map[f] = vector< FaceID >(); }
            face_steps_map[f].push_back( w.opp().face() );
            // prendo il prev, perché altrimenti parto una rotazione avanti
            vs.push_back( other.pos( w.prev().vertex( )));
            w = w.next();
        }
        FaceID  new_f               = m.add_face( vs );
        assert(new_f != InvalidFaceID );
        old_to_new_faces[f] = new_f;
    }
    // for each face of source - merge the faces using as reference the number of steps
    // needed to the walker to get the right edge
    for( auto f : other.faces( ))
    {
        FaceID  new_f       = old_to_new_faces[f];
        auto    neighbors   = face_steps_map[f];
        Walker  fw          = m.walker( new_f );
        for( int steps = 0; steps < neighbors.size(); ++steps )
        {
            //                    if( fw.opp().face() != InvalidFaceID ) continue;
            // find the neighor face
            FaceID old_f_neighbor = neighbors[steps];
            FaceID new_neighbor   = old_to_new_faces[old_f_neighbor];
            // find N that is the index at face_steps_map[old_f_neighbor that gives f
            int i = 0;
            for(  ; i < neighbors.size() && face_steps_map[old_f_neighbor][i] != f; ++i ); // loops without doing anything
            assert( face_steps_map[old_f_neighbor][i] == f );
            
            Walker neighbor_w = m.walker(new_neighbor);
            for( int j = 0; j < i; ++j ) { neighbor_w = neighbor_w.next(); }
            m.stitch_boundary_edges( neighbor_w.opp().halfedge(), fw.opp().halfedge( ));
            fw = fw.next();
        }
    }
    
    m.cleanup( remap );
    for( auto p : remap.vmap )
    {
        if( m_IDs_set.count( p.first )) { m_IDs_in_result.insert( p.second );   }
        else                            { other_IDs_in_result.insert( p.second ); }
    }
    
}



/// adds other to m, and gives in output the set of the IDs of the two original meshes inside the resulting mesh m
/// order is preserved
void add_manifold( Manifold &m, Manifold &other, set<VertexID> &m_IDs_in_result, set<VertexID> & other_IDs_in_result )
{
    IDRemap _;
    add_manifold(m, other, m_IDs_in_result, other_IDs_in_result, _ );
}
        
void test_delete ( Manifold &m, set<VertexID> &to_delete )
{
    for( auto v : to_delete )
    {
        m.remove_vertex(v);
    }
}

}}