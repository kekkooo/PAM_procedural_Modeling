//
//  structural_opeations.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "structural_operations.h"

#include <GEL/CGLA/Vec3d.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Operations/geometric_operations.h>
#include "polarize.h"

using namespace HMesh;
using namespace std;
using namespace CGLA;
using namespace Procedural::Geometry;
using namespace Procedural::Operations::Geometric;


namespace Procedural{
    namespace Operations{
        namespace Structural{

void add_branch ( HMesh::Manifold& m, HMesh::VertexID vid, int size, HMesh::VertexAttributeVector<int> &ring )
{
//    HMesh::VertexAttributeVector<int> ring( m.no_vertices(), 0 );
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
    // seleziona solo i punti   che non sono poli.
    //                          che non hanno vicini poli
    //                          che non hanno vicini selezionati
        
    map< HMesh::VertexID, CGLA::Vec3d > vert_pos;
    if( !is_pole(m, vid) )
    {
        Walker w  = m.walker( vid );
        ring[vid] = 1;
        
        for ( ; !w.full_circle(); w=w.circulate_vertex_ccw() )
        {
            if ( edge_info[w.halfedge()].is_rib() )
            {
                        ring[w.vertex()]    = 1;
                Walker  side_w              = w;
                for (int sides = 1; sides < size; ++sides )
                {
                    side_w = side_w.next().opp().next();
                    ring[side_w.vertex()] = 1;
                }
            }
            vert_pos[w.vertex()] = m.pos( w.vertex( ));
        }
        Vec3d centroid(0);
        polar_add_branch( m, ring );
    }
}
            
void cut_branch( Manifold& m, HalfEdgeID h, VertexID pole )
{
    //  0) find the ring of vertices
    vector< VertexID > vertices;
    
    // 1) save a vector from one of the pole's neighbors to the pole ( or I need all of them? )
    // need to save anche which of the vertices of the slitted ring has to be used as reference
    Vec3d       barycenter      = ring_barycenter( m, h, vertices );
    Vec3d       pole_dir        = m.pos(pole) - barycenter;
    // 2) slit that ring of vertices
    VertexAttributeVector<int> selection(m.no_vertices(), 0);
    for( VertexID v : vertices ) { selection[v] = 1; }
    HalfEdgeID hole_edge = m.slit_edges(selection);
    // 3) delete all the faces of the slitted component
    // this should cut also the sub_branches!!!!!!!!!!!! ( NOT SO EASY! )
    vector< VertexID > v_to_delete;
    v_to_delete.push_back(pole);
    Walker  w       = m.walker(pole);
    bool    done    = false;
    w       = m.walker( w.next().halfedge( ));
    
    while( !done )
    {
        // this just works for a single branch with no sub-branches
        /* starting from the pole
         for each outgoint branch
         mark the vertex to be deleted and move the walker to next.opp.next
         // each walker stops when reaches a vertex inside "vertices" ( put them in a SET! )
         */
        
        
        
        
        //        done = ( w.opp().face() == InvalidFaceID);
        //        for (; !w.full_circle(); w = w.next().opp().next()) {
        //            v_to_delete.push_back(w.vertex());
        //        }
        //        if (!done)
        //        {
        //            w = m.walker( w.opp().next().next().halfedge( ));
        //        }
    }
    for ( auto vtd : v_to_delete )
    {
        cout << "deleting vertex " << vtd << endl;
        m.remove_vertex(vtd);
    }
    // 4) rebuild the pole
    FaceID      f           = m.close_hole( hole_edge );
    VertexID    new_pole    = m.split_face_by_vertex( f );
    Vec3d       f_centroid  = ring_barycenter( m, hole_edge );
    m.pos( new_pole )       = f_centroid + pole_dir;
}

// this works only for cutting branches that are not part of a loop
// NEED TO KNOW how to understand if the ring is part of a loop
// ALSO need to understend in which direction should the walker move in order to reach the pole
void cut_branch ( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    // if not loop OR not between two joints
    //              find the right pole
    //              call cut_branch( m, h, pole )
    // else cut creating two poles with a simple policy
}
   
            
}}}
