//
//  graph_match.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 10/10/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "graph_match.h"
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>

//#include "geometric_properties.h"
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
        

void fill_graph( const PoleList &poles, const PoleInfoMap& poleInfo, GraphStruct &g, ManifoldToGraph &mtg ){
    // calculate costs save them into the graph
    for( GraphEdge e : g.arcs )
    {
        VertexID v1 = mtg.getVertexId( e.first ),
                 v2 = mtg.getVertexId( e.second );
        
        assert( poleInfo.count( v1 ) > 0 );
        assert( poleInfo.count( v2 ) > 0 );
        
        CGLA::Vec3d n1 = poleInfo.at( v1 ).geometry.normal,
                    n2 = poleInfo.at( v2 ).geometry.normal;
        n1.normalize();
        n2.normalize();
#warning consider using squared length
        double distance = ( poleInfo.at( v1 ).geometry.pos - poleInfo.at( v2 ).geometry.pos).length();
        double cos      = dot( n1, n2 ) + 1.0;
        // since cos is in the range [-1, 1] and I use differences, I apply an offset to the range [ 0, 2 ]
        
        assert( distance > 0 );
        assert( !isnan( cos ));
        g.setCost( e, make_pair( distance, cos ));
    }
}
        
        
void normalize_costs( GraphStruct &g1, GraphStruct g2, double max_distance ){
    // normalize edge cost values
    for( GraphEdge e : g1.arcs )
    {
        g1.setCost( e, std::make_pair( g1.getCost( e ).first / max_distance, g1.getCost( e ).second / 2.0 ));
#warning here I need a well thoughed assert
    }
    for( GraphEdge e : g2.arcs )
    {
        g2.setCost( e, std::make_pair( g2.getCost( e ).first / max_distance, g2.getCost( e ).second / 2.0 ));
#warning here I need a well thoughed assert
    }
}

void get_subsets( const MainStructure& main, const Module& module,
                  const std::vector< Match >& proposed,
                  std::vector< SubsetResult >& result, EdgeCost treshold ){
    
    size_t no_nodes = module.poleList.size();
    EdgeCost cost_sum = make_pair( 0.0, 0.0 );
    
    PoleList main_poles, module_poles;
    
    for( const auto& pole_and_vertex : proposed )
    {
        module_poles.push_back(  pole_and_vertex.first);
        main_poles.push_back( pole_and_vertex.second );
    }
    
    GraphStruct gm( no_nodes );
    GraphStruct gh( no_nodes );
    GraphStruct g( 0 );
    ManifoldToGraph mtg_module( module_poles );
    ManifoldToGraph mtg_main( main_poles );

    // fill graph costs for module
    fill_graph( module_poles, module.getPoleInfoMap(), gm, mtg_module );
    // fill graph costs for main
    fill_graph( main_poles, main.getPoleInfoMap(), gh, mtg_main );
    // build difference graph
    graphStruct_difference( gm, gh, g );
    
    SubsetResult s0;
    getSubsetResult( g, mtg_main, mtg_module, s0 );
    
    bool go_on = s0.cost < treshold;
    
    // iterate removing until the target is found
    for( size_t i = 0; i < no_nodes && go_on; ++i )
    {
        result.push_back( s0 );
        
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
        
        EdgeCost    max_cost = ( *total_cost.begin( )).second;
        GraphNode   choosen  = ( *total_cost.begin( )).first;
        
        // scegliere il vertice con il costo maggiore della star
        for( const auto& star_cost : total_cost )
        {
            if( star_cost.second > max_cost )
            {
                max_cost = star_cost.second;
                choosen  = star_cost.first;
            }
        }
        // rimuoverlo dal grafo
        g.RemoveNode( choosen );
        getSubsetResult( g, mtg_main, mtg_module, s0 );
    }
    
    assert( g.no_nodes() > 0 );
}
    
        
void getSubsetResult( const GraphStruct& g, const ManifoldToGraph& mtg_main, const ManifoldToGraph& mtg_module,
                      SubsetResult& result ){
    result.cost = make_pair( 0.0, 0.0 );
    result.matches.clear();
    for( GraphNode index : g.nodes )
    {
        if( g.exists( index )) {
            result.cost = result.cost + getStarTotalCost( g, index );
            result.matches.push_back( make_pair( mtg_module.getVertexId( index ), mtg_main.getVertexId( index )));
        }
    }
}

EdgeCost getGraphTotalCost( const GraphStruct& g ){
    EdgeCost cost_sum = make_pair( 0.0, 0.0 );
    // save the matches
    for( GraphNode index : g.nodes )
    {
        if( g.exists( index )) {
            cost_sum = cost_sum + getStarTotalCost( g, index );
        }
    }
    return cost_sum;
}
        
EdgeCost getStarTotalCost( const GraphStruct &g, GraphNode n )
{
    vector< GraphEdge > star;
    EdgeCost cost = make_pair( 0.0, 0.0 );
    // prendere la sua star
    g.getStar( n, star );
    // calcolare la somma dei costi della star
    for( GraphEdge e : star )
    {
        cost = cost + g.getCost( e );
    }
    return cost;
}
        
        
GraphNode remove_most_expensive_node( GraphStruct &g )
{
    map< GraphNode, EdgeCost > total_cost;
    // per ogni vertice
    for( size_t idx : g.nodes )
    {
        if( g.exists( idx ))
        {
            total_cost[idx] = getStarTotalCost( g, idx );
        }
    }
    
    EdgeCost    max_cost = (*total_cost.begin()).second;
    GraphNode   choosen  = (*total_cost.begin()).first;
    
    // scegliere il vertice con il costo maggiore della star
    for( const auto& star_cost : total_cost )
    {
        if( star_cost.second > max_cost )
        {
            max_cost = star_cost.second;
            choosen  = star_cost.first;
        }
    }
    // rimuoverlo dal grafo
    g.RemoveNode( choosen );
    return choosen;
}
        
/********************************************/
/*                  UTILITIES               */
/********************************************/
std::ostream& graph_print (std::ostream &out, EdgeCost cost)
{
    out << "(" << cost.first << ", " << cost.second << ")";
    return out;
}
   
        
std::ostream& graph_print (std::ostream &out, GraphNode &node)
{
    out << "[" << (long)node << "]";
    return out;
}

        
std::ostream& graph_print (std::ostream &out, GraphEdge &edge)
{
    out << "[" << (long)edge.first << ", " << (long)edge.second << "]";
    return out;
}
        
        
std::ostream& operator<< (std::ostream &out, GraphStruct &g)
{
        out << " the graph has : " << g.no_nodes() << " nodes and " << g.arcs.size() << " edges " << std::endl;
        out << "the nodes are :" << std::endl;
        for( GraphNode n : g.nodes )
        {
            if(g.exists(n))
            {
                graph_print( out, n )  << std::endl;
            }
        }
        out << "the arcs are :" << std::endl;
        for( const auto& e : g.arcs )
        {
            graph_print(out, e) ;
            out << " ==> ";
            graph_print(out, g.getCost(e))  << std::endl ;
//            out <<  std::endl;
        }
        out << "================END=====================" << std::endl << std::endl;
    //
    return out;
}


}}