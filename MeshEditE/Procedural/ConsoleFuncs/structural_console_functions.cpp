//
//  structural_console_functions.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "structural_console_functions.h"
#include <GEL/GLGraphics/MeshEditor.h>
//#include <strstream>
//#include <istream>
//#include <fstream>
//#include <string>
//#include <regex>
//#include <set>
//#include <queue>
//#include "Test.h"
#include "polarize.h"
#include <MeshEditE/Procedural/Operations/structural_operations.h>

using namespace GLGraphics;
using namespace std;
using namespace HMesh;
using namespace Procedural::Operations::Structural;

// THE LOGIC OF THIS THING SHOULD BE MOVED INTO ANOTHER CLASS
void console_test_add_branch( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m           = me->active_mesh();
    int         size        = 1;
    vector< VertexID > selected;

    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> size;
    }
    
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit )
    {
        if( me->get_vertex_selection()[*vit] ) selected.push_back(*vit);
    }
    // seleziona solo i punti   che non sono poli.
    //                          che non hanno vicini poli
    //                          che non hanno vicini selezionati
    for( VertexID vid : selected )
    {
        if( me->get_vertex_selection()[vid] )
            { add_branch( m, vid, size, me->get_vertex_selection( )); }
    }
    
    me->active_visobj().clear_selection();
    m.cleanup();
    
    cout << "add branches to selected vertices" << endl;
}

void console_test_cut_branch( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    bool        done            = false;
    Manifold&   m               = me->active_mesh();
    VertexID    selected_pole   = InvalidVertexID;
    HalfEdgeID  selected_he     = InvalidHalfEdgeID;
    
    vector< VertexID >                      selected;
    typedef vector< VertexID >::iterator    vertexID_iter;
    HalfEdgeAttributeVector<EdgeInfo>       edge_info = label_PAM_edges( m );

    // consider only selected vertices that are not poles
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit] )
        {
            selected.push_back(*vit);
        }
    }
    
    // take a couple of selected vertices that define a spine edge
    //    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end() && !done; ++vit)
    for( vertexID_iter vit = selected.begin();
        vit != selected.end() && ( !done || selected_pole == InvalidVertexID ); ++vit)
    {
        
        if( is_pole(m, *vit ))
        {
            selected_pole = *vit;
        }
        else if( !done )
        {
            Walker w = m.walker( *vit );
            for (; !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                assert(*vit != w.vertex());
                // no need to check if it is not a pole because from a pole you can't have ribs
                if( me->get_vertex_selection()[w.vertex()] && edge_info[w.halfedge()].is_rib()  )
                {
                    selected_he = w.halfedge();
                    done = true;
                }
            }
        }
    }
    assert(selected_he != InvalidHalfEdgeID);
    assert( selected_pole != InvalidVertexID);
    //    split_ring_of_quads(m, selected_he);
    cut_branch(m, selected_he, selected_pole);
    m.cleanup();
}


namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_structural_console_funcs(GLGraphics::MeshEditor* me)
        {
            me->register_console_function( "test.structure.add_branch", console_test_add_branch, "test.structure.add_branch" );
            
            me->register_console_function( "test.structure.cut_branch", console_test_cut_branch, "test.structure.cut_branch" );
        }
}}