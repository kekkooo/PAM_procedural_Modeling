//
//  graph_match.h
//  MeshEditE
//
//  Created by Francesco Usai on 10/10/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__graph_match__
#define __MeshEditE__graph_match__

#include <stdio.h>
#include <stack>
#include <ostream>
#include <GEL/HMesh/Manifold.h>

#define EDGE_COST_EPS 0.00000001

namespace Procedural{
    namespace GraphMatch{
        
typedef size_t                     GraphNode;
typedef std::pair<size_t, size_t>  GraphEdge;
typedef std::pair<double, double>  EdgeCost;

// should define here the operators +,-, >, <, == for the EdgeCost
inline EdgeCost operator +( const EdgeCost& l, const EdgeCost& r )
        { return std::make_pair( l.first + r.first, l.second + r.second ); };

inline EdgeCost operator -( const EdgeCost& l, const EdgeCost& r )
        { return std::make_pair( fabs( l.first - r.first ), fabs( l.second - r.second )); };
        
inline bool is_equal( double l, double r ) { return fabs( l - r ) < EDGE_COST_EPS; }
inline bool operator ==( const EdgeCost& l, EdgeCost& r ) { return is_equal( l.first, r.first ) && is_equal( l.second, r.second ); }
inline bool operator !=( const EdgeCost& l, EdgeCost& r ) { return !( l == r ); }
        
inline bool operator <( const EdgeCost& l, const EdgeCost& r)
{
    if( l.first < r.first && l.second < r.second )
    {
        return true;
    }
    else if( is_equal( l.first, r.first ))
    {
        return ( l.second < r.second );
    }
    else if( is_equal( l.second, r.second ))
    {
        return ( l.first < r.first );
    }
    else { return false; }
}
        
inline bool operator <=( const EdgeCost& l, const EdgeCost& r)
{
    return (( l < r ) || ( l == r ));
}
        
inline bool operator >( const EdgeCost& l, const EdgeCost& r)
{
    return !( l <= r );
}

inline bool operator >=( const EdgeCost& l, const EdgeCost& r)
{
    return !( l < r );
}
        
GraphEdge build_edge( GraphNode n1, GraphNode n2 );
        
struct ManifoldToGraph
{
private:
    std::map< HMesh::VertexID, GraphNode > id_to_node;
    std::map< GraphNode, HMesh::VertexID > node_to_id;
public:
    inline GraphNode        getNode( HMesh::VertexID v )   { assert( id_to_node.count(v) > 0 ); return id_to_node[v]; }
    inline HMesh::VertexID  getId  ( GraphNode n )         { assert( node_to_id.count(n) > 0 ); return node_to_id[n]; }
    ManifoldToGraph( const std::vector<HMesh::VertexID> &vs )
    {
        GraphNode curr = 0;
        for( HMesh::VertexID v : vs )
        {
            id_to_node[v]       = curr;
            node_to_id[curr]    = v;
            ++curr;
        }
    }
};

/// super simple graph structure, manages only complete graphs
/// nodes are numbered from 0 to no_nodes - 1
/// initialization is super-simple as well, just feed the number of nodes
/// it will create the complete graph, then set the various costs.
/// if a cost is not set it will be the maximum
/// nodes cannot be added to this structure
/// nodes can be deleted from this structure, simply they will be set as -1 at the corresponding index
struct GraphStruct
{
private :
    std::map< GraphEdge, EdgeCost > costs;
    size_t                          max_node_idx;
    
public :
    std::vector<GraphNode> nodes;
    std::vector<GraphEdge> arcs;
    
    inline bool exists( GraphNode n ) { return ( n <= max_node_idx && nodes[n] == n ); }
    inline bool exists( GraphEdge e ) { return ( exists( e.first ) && exists( e.second )); }
    
    inline size_t no_nodes()
    {
        size_t count = 0;
        for ( auto n : nodes ) { if( exists(n)) ++count; }
        return count;
    }
    
    inline void setCost( GraphEdge e, EdgeCost cost )
    {
        assert( exists( e ));
        costs[e] = cost;
    }
    
    inline void setCost( GraphNode n1, GraphNode n2, EdgeCost cost )
    {
        assert( exists( n1 ));
        assert( exists( n2 ));
        costs[build_edge( n1, n2 )] = cost;
    }
    
    inline EdgeCost getCost( GraphEdge e )
    {
        assert( exists( e ));
        return costs[e];
    }
    
    inline EdgeCost getCost( GraphNode n1, GraphNode n2 )
    {
        assert( exists( n1 ));
        assert( exists( n2 ));
        return costs[build_edge( n1, n2 )];
    }
    
    inline void getStar( GraphNode n, std::vector< GraphEdge > &star )
    {
        assert( exists( n ));
        for( size_t i = 0; i < max_node_idx; ++i )
        {
            if( exists( i ))
            {
                if( i < n ) { star.push_back( std::make_pair( i, n )); }
                if( i > n ) { star.push_back( std::make_pair( n, i )); }
            }
        }
    }
    
    GraphStruct(){}
    
    GraphStruct( size_t no_nodes )
    {
        assert( no_nodes >= 0 );
        max_node_idx = no_nodes - 1;
        double ub   = std::numeric_limits<double>::max();
        EdgeCost c  = std::make_pair( ub, ub );
        // add nodes
        for( size_t i = 0; i < no_nodes; ++i )
        {
            nodes.push_back( (int)i );
        }
        // add arcs
        for( size_t i = 0; i <= no_nodes; ++i )
        {
            for( size_t j = i+1; j < no_nodes; ++j )
            {
                GraphEdge e = build_edge( (int)i, (int)j );
                arcs.push_back( e );
                costs[e] = c;
            }
        }
    }
    
    /// removing a node of the graph involves :
    /// putting -1 in its corresponding position in the node's vector
    /// removing all of its arcs
    /// putting 0 as cost for each of its arcs
    void RemoveNode( size_t node )
    {
        nodes[node] = -1;
        assert( !exists( node ));
        // save all the arcs where
        typedef std::vector<GraphEdge>::iterator    arc_iter;
        std::stack<arc_iter>                        to_delete;
        // mark every arc outgoing from the node to be deleted
        for( arc_iter it = arcs.begin(); it != arcs.end(); ++it )
        {
            if( it->first == node || it->second == node )
                to_delete.push(it);
        }
        // remove them ( popping from the stack will guarantee of deleting them from back to front
        while( !to_delete.empty( ))
        {
            costs[(*to_delete.top())] = std::make_pair( 0.0, 0.0 );
            arcs.erase( to_delete.top( ));
            to_delete.pop();
        }
    
    }
};


std::ostream& graph_print (std::ostream &out, EdgeCost cost);
std::ostream& graph_print (std::ostream &out, GraphNode &node);
std::ostream& graph_print (std::ostream &out, GraphEdge &edge);
std::ostream& operator<< (std::ostream &out, GraphStruct &g);



void        graphStruct_difference(         GraphStruct &g1, GraphStruct &g2, GraphStruct &g );
GraphNode   remove_most_expensive_node (    GraphStruct &g );
EdgeCost    getStarTotalCost(               GraphStruct &g, GraphNode n );

EdgeCost get_best_subset( HMesh::Manifold &host,  std::vector<HMesh::VertexID> &aps,
                          HMesh::Manifold &module, std::vector<HMesh::VertexID> &poles,
                          std::vector< size_t > &selected_indices, size_t target );

void     fill_graph    ( HMesh::Manifold &m, std::vector<HMesh::VertexID> &vs, GraphStruct &g );

}}
#endif /* defined(__MeshEditE__graph_match__) */
