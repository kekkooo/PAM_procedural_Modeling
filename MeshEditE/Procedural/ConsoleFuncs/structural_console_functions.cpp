//
//  structural_console_functions.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "structural_console_functions.h"
#include <GEL/GLGraphics/MeshEditor.h>
#include "polarize.h"
#include <MeshEditE/Procedural/Operations/structural_operations.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MeshEditE/Procedural/Helpers/module_alignment.h>
#include <GEL/HMesh/obj_load.h>


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

void console_test_glue_poles( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m               = me->active_mesh();
    
    typedef vector< VertexID >::iterator    vertexID_iter;
//    HalfEdgeAttributeVector<EdgeInfo>       edge_info = label_PAM_edges( m );
    vector< VertexID > poles;
    
    // consider only selected vertices that are not poles
    for( VertexIDIterator vit = m.vertices_begin(); vit != m.vertices_end() && poles.size() < 2; ++vit)
    {
        if (me->get_vertex_selection()[*vit] )
        {
            poles.push_back(*vit);
        }
    }
    if( poles.size() == 2 )
        glue_poles( m, poles[0], poles[1] );
    
    m.cleanup();
}

void console_test_add_branch_on_high_angles( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    Manifold&   m       = me->active_mesh();
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges( m );
    VertexAttributeVector<double> angles;
    
    Procedural::Structure::LabelJunctions( m, edge_info );
    Procedural::Geometry::dihedral_angles( m, edge_info, angles );

    for( auto vid : m.vertices( ))
    {
        cout << vid << "has cosine " << angles[vid] << endl;
        if (angles[vid] > 0.0 )
        {
            me->active_visobj().clear_selection();
            me->get_vertex_selection()[vid] = 1;
            add_branch( m, vid, 1, me->get_vertex_selection( ));
        }
    }
    me->active_visobj().clear_selection();
    m.cleanup();
}

void console_test_add_module( MeshEditor *me, const std::vector< std::string > &args )
{
    me->save_active_mesh();
    stringstream    oss;
    Manifold        module;
    Manifold&       host                = me->active_mesh();
    string          filename            = "test_match.obj";
    int             target_glueings     = 2;
    std::vector<Procedural::Match> matches;

    if( args.size() > 0 ){
        istringstream a0( args[0] );
        a0 >> filename;
    }
    
    if( args.size() > 1 ){
        istringstream a1( args[1] );
        a1 >> target_glueings;
    }
    
    oss << "/Users/francescousai/Documents/Dottorato/Visiting/Results/8v_poles/" << filename;
    obj_load( oss.str(), module );
//    Procedural::Helpers::ModuleAlignment::AddModule( host, module, target_glueings, matches );
    for( Procedural::Match m : matches )
    {
        me->get_vertex_selection()[m.first]     = 1;
        me->get_vertex_selection()[m.second]    = 1;
    }
    
}

namespace Procedural{
    namespace ConsoleFuncs{
        
        void register_structural_console_funcs(GLGraphics::MeshEditor* me)
        {
            me->register_console_function( "test.structure.add_branch", console_test_add_branch, "test.structure.add_branch" );
            
            me->register_console_function( "test.structure.cut_branch", console_test_cut_branch, "test.structure.cut_branch" );

            me->register_console_function( "test.structure.remove_branch", console_test_remove_branch, "test.structure.remove_branch" );
            
            me->register_console_function( "test.structure.glue_poles", console_test_glue_poles, "test.structure.glue_poles" );
            
            me->register_console_function( "test.structure.add_branch_on_high_angles", console_test_add_branch_on_high_angles, "test.structure.add_branch_on_high_angles" );

            me->register_console_function( "test.add_module", console_test_add_module, "console_test_add_module" );
            
            
            
        }
}}