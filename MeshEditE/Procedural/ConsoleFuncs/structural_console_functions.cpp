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
    Manifold&   m               = me->active_mesh();
    
    typedef vector< VertexID >::iterator    vertexID_iter;
    HalfEdgeAttributeVector<EdgeInfo>       edge_info = label_PAM_edges( m );

    // consider only selected vertices that are not poles
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit] )
        {
            cut_branch(m, *vit, edge_info);
        }
    }
    m.cleanup();
}

void console_test_remove_branch( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m               = me->active_mesh();
    
    typedef vector< VertexID >::iterator    vertexID_iter;
    HalfEdgeAttributeVector<EdgeInfo>       edge_info = label_PAM_edges( m );
    
    // consider only selected vertices that are not poles
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end(); ++vit)
    {
        if (me->get_vertex_selection()[*vit] )
        {
            remove_branch(m, *vit, edge_info);
        }
    }
    m.cleanup();
}

namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_structural_console_funcs(GLGraphics::MeshEditor* me)
        {
            me->register_console_function( "test.structure.add_branch", console_test_add_branch, "test.structure.add_branch" );
            
            me->register_console_function( "test.structure.cut_branch", console_test_cut_branch, "test.structure.cut_branch" );

            me->register_console_function( "test.structure.remove_branch", console_test_remove_branch, "test.structure.remove_branch" );
            
        }
}}