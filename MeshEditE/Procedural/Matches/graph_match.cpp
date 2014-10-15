//
//  graph_match.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 10/10/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "graph_match.h"
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>

using namespace std;
using namespace HMesh;

namespace Procedural{
    namespace GraphMatch{
        
void graphStruct_difference( GraphStruct &g1, GraphStruct &g2, GraphStruct &g )
{
    // they should have the same node set and arc set
    // will check only the number of nodes ( assumptions can be made accordingly to the constructor )
    assert( g1.nodes.size() == g1.nodes.size() );
    size_t no_nodes = g1.nodes.size();
    g = GraphStruct( no_nodes );
    
    for( GraphEdge e : g.arcs )
    {
        g.setCost( e, g1.getCost( e ) - g2.getCost( e ));
    }
}

GraphEdge build_edge( GraphNode n1, GraphNode n2 )
{
    assert( n1 != n2 );
    if( n1 < n2 ) return std::make_pair( n1, n2 );
    else          return std::make_pair( n2, n1 );
    
}
        
void fill_graph ( Manifold &m, vector<HMesh::VertexID> &vs, GraphStruct &g, ManifoldToGraph &mtg )
{
    vector< CGLA::Vec3d > pos, normals;
    for( VertexID v : vs )
    {
        pos.push_back( m.pos( v ));
        normals.push_back( Procedural::Geometry::vertex_normal( m, v ));
    }
    // calculate costs save them into the graph
    for( GraphEdge e : g.arcs )
    {
        VertexID v1 = mtg.getId( e.first ),
                 v2 = mtg.getId( e.second );
        CGLA::Vec3d n1 = Procedural::Geometry::vertex_normal( m, v1 ),
                    n2 = Procedural::Geometry::vertex_normal( m, v2 );

        double distance = ( m.pos( v1 ) - m.pos( v2 )).length();
        double angle    = acos( CGLA::dot( n1, n2 ));
        g.setCost( e, make_pair( distance, angle ));
    }
}
    
        
EdgeCost get_best_subset( Manifold &host,   vector< VertexID > &aps,
                          Manifold &module, vector< VertexID > &poles,
                          vector< size_t > &selected_indices, size_t target )
{
    assert( aps.size() == poles.size( ));
    assert( target <= aps.size());
    size_t no_nodes = aps.size();
    EdgeCost cost_sum = make_pair( 0.0, 0.0 );
    GraphStruct gm( no_nodes );
    GraphStruct gh( no_nodes );
    GraphStruct g( 0 );
    ManifoldToGraph mtg_module( poles );
    ManifoldToGraph mtg_host( aps );

    // fill graph costs for module
    fill_graph( module, poles, gm, mtg_module );
    // fill graph costs for host
    fill_graph( host, aps, gh, mtg_host );
    // build difference graph
    graphStruct_difference( gm, gh, g );
    // iterate removing until the target is found
    size_t iterations = no_nodes - target;
    for( size_t i = 0; i < iterations; ++i )
    {
        map< GraphNode, EdgeCost > total_cost;
        // per ogni vertice
        for( size_t idx : g.nodes )
        {
            if( g.exists( idx ))
            {
                vector< GraphEdge > star;
                // prendere la sua star
                g.getStar( idx, star );
                total_cost[idx] = make_pair( 0.0, 0.0 );
                // calcolare la somma dei costi della star
                for( GraphEdge e : star )
                {
                    total_cost[idx] = total_cost[idx] + g.getCost( e );
                }
                
            }
        }
        
        EdgeCost    max_cost = (*total_cost.begin()).second;
        GraphNode   choosen  = (*total_cost.begin()).first;
        
        // scegliere il vertice con il costo maggiore della star
        for( auto star_cost : total_cost )
        {
            if( star_cost.second > max_cost )
            {
                max_cost = star_cost.second;
                choosen  = star_cost.first;
            }
        }
        // rimuoverlo dal grafo
        g.RemoveNode( choosen );
    }
    return cost_sum;
}


}}