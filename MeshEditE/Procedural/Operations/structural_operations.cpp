//
//  structural_opeations.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "structural_operations.h"

#include "Algorithms.h"

#include <GEL/CGLA/Vec3d.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Operations/geometric_operations.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <MeshEditE/Procedural/Operations/Algorithms.h>
#include "polarize.h"

using namespace HMesh;
using namespace std;
using namespace CGLA;
using namespace Procedural::Geometry;
using namespace Procedural::Operations::Geometric;
using namespace Procedural::Structure;


namespace Procedural{
    namespace Operations{
        namespace Structural{

VertexID add_branch ( HMesh::Manifold& m, HMesh::VertexID vid, int size, HMesh::VertexAttributeVector<int> &ring )
{
//    HMesh::VertexAttributeVector<int> ring( m.no_vertices(), 0 );
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    
    // seleziona solo i punti   che non sono poli.
    //                          che non hanno vicini poli
    //                          che non hanno vicini selezionati
        
    map< HMesh::VertexID, CGLA::Vec3d > vert_pos;
    if( !is_pole( m, vid ) && !is_2_neighbor_of_pole( m,vid ))
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
        return polar_add_branch( m, ring );
    }
    return InvalidVertexID;
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
            
void cut_branch ( Manifold& m, VertexID v, HalfEdgeAttributeVector<EdgeInfo> edge_info )
{
    vector< VertexID > vertices;
    HalfEdgeID he = get_ring_vertices( m, v, edge_info, vertices, false );
    // early termination if the ring is not found
    if ( vertices.size() == 0 ) return;
    // save two references to the two sides ( from v )
    Walker w = m.walker(he);
    HalfEdgeID rib1 = w.next().next().halfedge();
    HalfEdgeID rib2 = w.opp().next().next().halfedge();
    for( VertexID vid : vertices ) { m.remove_vertex(vid); }
    assert( m.walker(rib1).face() == InvalidFaceID );
    assert( m.walker(rib2).face() == InvalidFaceID );
    FaceID f1 = m.close_hole(rib1);
    FaceID f2 = m.close_hole(rib2);
    m.split_face_by_vertex(f1);
    m.split_face_by_vertex(f2);
    
    // NEED TO DECIDE WHAT TO DO IN SOME CASES ( avoid separate components, etc )
}
            
void remove_branch ( HMesh::Manifold& m, HMesh::VertexID pole, HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info )
{
    // check if m.walker(pole).next().halfedge è di tipo junction
    if ( !is_pole( m, pole )) return;
    
    LabelJunctions( m, edge_info );
    
    if ( !edge_info[ m.walker(pole).next().opp().halfedge()].is_junction() ) return;

    typedef pair< HalfEdgeID, HalfEdgeID > ToStitchPair;
    vector< ToStitchPair > to_stitch;

    // save the pair of edges that must be stitched together
    Walker w = m.walker(pole);
    for( ; !w.full_circle() && !is_singularity(m, w.vertex()); w = w.circulate_vertex_ccw());
    w = w.next();

    assert( edge_info[w.halfedge()].is_junction( ));
    // remove pole
    m.remove_vertex( pole );
    Walker wback = w.prev();
        
    assert( is_singularity(m, wback.vertex( )));
    do
    {
        to_stitch.push_back( make_pair( w.halfedge(), wback.halfedge( )));
        w       = w.next();
        wback   = wback.prev();
        assert( w.face()     == InvalidFaceID );
        assert( wback.face() == InvalidFaceID );
        
    } while( valency( m, w.vertex()) != 5 );
//    to_stitch.push_back( make_pair( w.halfedge(), wback.halfedge( )));
    
    // stich pair of edges ( note that at the vertices with valence 6 you should stich
    // togheter the junction outgoing edge and its prev ( around the hole )
    for( ToStitchPair tsp : to_stitch )
    {
        // debug
        VertexID vfirst  = m.walker(tsp.first).vertex();
        VertexID vsecond = m.walker(tsp.second).vertex();
        
        
        m.stitch_boundary_edges( tsp.first, tsp.second );
        cout << " is first? # "
             << m.in_use(vfirst) << " or is second? # "
             << m.in_use(vsecond) << endl;;
        
        //debug
        
    }
    
    // trovare un modo per risistemare le posizioni dei vertici
    Procedural::Operations::Algorithms::along_spines(m, edge_info);
}
            
            
void glue_poles ( Manifold& m, VertexID pole1, VertexID pole2 )
{
    // early termination in case the input vertices are not poles
    if( !is_pole( m, pole1 )) return;
    if( !is_pole( m, pole2 )) return;
    // get valency of the input vertices
    int         pole1_valency  = valency( m, pole1 ),
                pole2_valency  = valency( m, pole2 );
    // SKIP THE FOLLOWING IF THE VALENCIES ARE EQUAL
    // this part subdivides the low valency pole  in order to match valencies
    if( pole1_valency != pole2_valency )
    {
        // choose the one that has lowest valency ( say low_val_pole )
        VertexID    low_val_pole   = pole1_valency < pole2_valency ? pole1 : pole2;
        VertexID    hi_val_pole    = low_val_pole == pole1         ? pole2 : pole1;
        int         lvpole_valency = low_val_pole == pole1         ? pole1_valency : pole2_valency;
        int         hvpole_valency = low_val_pole == pole1         ? pole2_valency : pole1_valency;
        int         diff           = hvpole_valency - lvpole_valency;
        assert( diff > 0 );
        // subdivide it for a number of times equal to the difference in valence
        HalfEdgeID starter = m.walker( low_val_pole ).next().halfedge();
        // could happen that the difference in valence is greater than the low_val_pole valency
        bool done = false;

        do
        {
            vector< VertexID >      vs;
            vector< HalfEdgeID >    hes;
            ring_vertices_and_halfedges( m, starter, vs, hes );
            
            size_t limit = diff > hes.size() ? hes.size() : diff;
            diff -= limit;
            size_t pace = hes.size() / limit;
            done = ( diff <= 0 );
            
            for( int i = 0; i < limit; i++ )
            {
                size_t index = i * pace;
                assert(i < hes.size());
                split_from_pole_to_pole(m, hes[index]);
            }
            
        } while( !done );
    }
    // get the two halfedge sequences
    vector< VertexID >      vs_pole1;
    vector< HalfEdgeID >    hes_pole1;
    ring_vertices_and_halfedges( m, m.walker(pole1).next().halfedge(), vs_pole1, hes_pole1 );

    vector< VertexID >      vs_pole2;
    vector< HalfEdgeID >    hes_pole2;
    ring_vertices_and_halfedges( m, m.walker(pole2).next().halfedge(), vs_pole1, hes_pole2 );
    assert(hes_pole1.size() == hes_pole2.size());

    // find adequate matchings between edges
    // use as matching policy the minimum angle between the vectors defined by :
    // - the two poles
    // - the two to_vertex of two compared halfedges
    Vec3d       pole_pole               = m.pos(pole1) - m.pos(pole2);
    double      min_angle_sum           = numeric_limits<double>::max();
    double      min_angle_diff          = numeric_limits<double>::max();
    HalfEdgeID  candidate_from_pole1    = InvalidHalfEdgeID;
    HalfEdgeID  candidate_from_pole2    = InvalidHalfEdgeID;
    int         candidate_id_from_pole1 = -1;
    int         candidate_id_from_pole2 = -1;
    
//    reverse(hes_pole2.begin(), hes_pole2.end());
    
    //for( HalfEdgeID p1_he : hes_pole1 )
    for( int i = 0; i < hes_pole1.size(); ++i )
    {
        HalfEdgeID p1_he = hes_pole1[i];
        Vec3d p1_a   = m.pos( m.walker( p1_he ).vertex());
        Vec3d p1_b = m.pos( m.walker( p1_he ).prev().vertex());
        
        for( int j = 0; j < hes_pole2.size(); ++j )
        {
            HalfEdgeID  p2_he  = hes_pole2[j];
            Vec3d       p2_b   = m.pos( m.walker( p2_he ).vertex());
            Vec3d       p2_a   = m.pos( m.walker( p2_he ).prev().vertex());
            // calculate the vectors
            Vec3d       a_dir       = p1_a - p2_a;
            Vec3d       b_dir       = p1_b - p2_b;
            double      angle_a     = angle( a_dir, pole_pole );
            double      angle_b     = angle( b_dir, pole_pole );
            double      angle_sum   = angle_a + angle_b;
            double      angle_diff  = fabs( angle_a - angle_b );
            
            cout << "current min " << min_angle_sum << " # " << min_angle_diff  << endl;
            cout << "latest calc " << angle_sum     << " # " << angle_diff      << endl;

            if( angle_sum <= min_angle_sum && angle_diff <= min_angle_diff )
            {
                min_angle_sum           = angle_sum;
                min_angle_diff          = angle_diff;
                candidate_from_pole1    = p1_he;
                candidate_from_pole2    = p2_he;
                candidate_id_from_pole1 = i;
                candidate_id_from_pole2 = j;
            }
        }
    }
    int p1_offset = candidate_id_from_pole1;
    int p2_offset = candidate_id_from_pole2;
    
    cout << p1_offset << " # " << p2_offset << endl;
    
    // delete poles
    m.remove_vertex(pole1);
    m.remove_vertex(pole2);
    
    Walker pole1_chain_walker = m.walker(hes_pole1[p1_offset]);
    Walker pole2_chain_walker = m.walker(hes_pole2[p2_offset]);
    
    hes_pole1.clear();
    hes_pole2.clear();
    
    while(!pole1_chain_walker.full_circle())
    {
        hes_pole1.push_back( pole1_chain_walker.halfedge( ));
        hes_pole2.push_back( pole2_chain_walker.halfedge( ));
        pole1_chain_walker = pole1_chain_walker.next();
        pole2_chain_walker = pole2_chain_walker.prev();
    }
    assert( pole2_chain_walker.full_circle( ));
    assert( hes_pole1.size() == hes_pole2.size( ));
    
    vector< pair< HalfEdgeID, HalfEdgeID > > he_to_stitch;
    HalfEdgeID hook_for_splitting = InvalidHalfEdgeID;
    // save pairs of edges that have to be stitched
    for( int i = 0; i < hes_pole1.size(); i++ )
    {
        
        vector<Vec3d> new_vs;
        new_vs.push_back( m.pos( m.walker( hes_pole1[i] ).prev().vertex() )     );
        new_vs.push_back( m.pos( m.walker( hes_pole1[i] ).vertex() )            );
        new_vs.push_back( m.pos( m.walker( hes_pole2[i] ).prev().vertex() )     );
        new_vs.push_back( m.pos( m.walker( hes_pole2[i] ).vertex() )            );
        
        FaceID new_f = m.add_face( new_vs );
        
        Walker new_f_walker = m.walker( new_f );
        HalfEdgeID  h1 = new_f_walker.opp().halfedge(),
                    h3 = new_f_walker.next().next().opp().halfedge();
        if( hook_for_splitting == InvalidHalfEdgeID ) hook_for_splitting = new_f_walker.next().halfedge();
        
        he_to_stitch.push_back( make_pair( hes_pole1[i], h1 ));
        he_to_stitch.push_back( make_pair( hes_pole2[i], h3 ));
        
    }
    
    assert(pole2_chain_walker.full_circle());

    for( int i = 0; i < he_to_stitch.size(); i++ )
    {
        auto hep = he_to_stitch[i];
        assert( m.stitch_boundary_edges( hep.first, hep.second ));
    }
    
    
    
    split_ring_of_quads( m, hook_for_splitting );
}
    
    
   
            
}}}
