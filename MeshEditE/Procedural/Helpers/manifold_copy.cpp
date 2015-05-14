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
        




/// adds other to m, and gives in output the set of the IDs of the two original meshes inside the resulting mesh m
/// order is preserved
void add_manifold( Manifold &m, Manifold &other, set<VertexID> &m_IDs_in_result, set<VertexID> & other_IDs_in_result )
{
    // save the set of the original vertices IDs of destination
    set<VertexID>                   m_IDs_set;
    map< FaceID, vector<FaceID> >   face_steps_map;
    map<FaceID, FaceID>             old_to_new_faces;
    IDRemap                         remap;
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
        else                            { other_IDs_in_result.insert(p.second); }
    }
}
        
void test_delete ( Manifold &m, set<VertexID> &to_delete )
{
    for( auto v : to_delete )
    {
        m.remove_vertex(v);
    }
}

}}