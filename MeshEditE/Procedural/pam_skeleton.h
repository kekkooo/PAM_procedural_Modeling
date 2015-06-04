//
//  pam_skeleton.h
//  MeshEditE
//
//  Created by Francesco Usai on 01/06/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef MeshEditE_pam_skeleton_h
#define MeshEditE_pam_skeleton_h

#include <map>
#include <vector>
#include <queue>
#include <set>
#include <fstream>
#include <iostream>

#include <GEL/HMesh/Manifold.h>

namespace Procedural{
    typedef size_t NodeID;
    typedef size_t BoneID;
    
    enum SkelNodeType{ SNT_Pole, SNT_Rib, SNT_Junction };
    
    struct SkelNode{
        NodeID                  ID;
        SkelNodeType            type    = SNT_Rib;
        CGLA::Vec3d             pos;
        double                  radius  = 0.0;
        std::set<NodeID>        neighbors;
        BoneID                  boneID;
    };
    
    struct SkelBone{
        BoneID                  ID;
        std::vector<NodeID>   nodes;
    };
    
    struct Skeleton{
    private:
        
        void openBone( NodeID id ){
            assert( id >= 0 && id < nodes.size() );
            assert( nodes[id].type == SNT_Pole || nodes[id].type == SNT_Junction );
            SkelBone b;
            b.ID = bones.size();
            b.nodes.push_back( id );
            bones.push_back( b );
        }
        
        void buildNeighborhood( BoneID id ){
            const SkelBone& b = bones[id];
            nodes[b.nodes.front()].boneID = b.ID;
            for( int i = 1; i < b.nodes.size() - 1; ++i ){
                NodeID prev = b.nodes[i-1],
                curr = b.nodes[i],
                next = b.nodes[i+1];
                
                nodes[curr].neighbors.insert( prev );
                nodes[curr].neighbors.insert( next );
                nodes[prev].neighbors.insert( curr );
                nodes[next].neighbors.insert( curr );
                nodes[curr].boneID = b.ID;
            }
            nodes[b.nodes.back()].boneID = b.ID;
        }
        
        void closeBone( NodeID id ){
            assert( bones.back().nodes.size() > 0 );
            SkelNodeType t_front = nodes[bones.back().nodes.front()].type;
            assert( t_front == SNT_Pole || t_front == SNT_Junction );
            assert( id >= 0 && id < nodes.size() );
            assert( nodes[id].type == SNT_Pole || nodes[id].type == SNT_Junction );
            
            bones.back().nodes.push_back( id );
            buildNeighborhood( bones.back().ID );
        }
        
        void addToCurrentBone( NodeID id ){
            assert( id >= 0 && id < nodes.size() );
            assert( nodes[id].type != SNT_Pole );
            assert( nodes[id].type != SNT_Junction );
            
            bones.back().nodes.push_back( id );
        }
        
        void addToBone( BoneID b_id, NodeID n_id ){
            // ids exist
            assert( b_id >= 0 && b_id < bones.size() );
            assert( n_id >= 0 && n_id < nodes.size() );
            // bone is not closed
            assert( nodes[bones[b_id].nodes.back()].type != SNT_Pole );
            assert( nodes[bones[b_id].nodes.back()].type != SNT_Junction );
            // node is not an endpoint
            assert( nodes[n_id].type != SNT_Pole );
            assert( nodes[n_id].type != SNT_Junction );

            bones[b_id].nodes.push_back( n_id );
        }
        
        void closeBone( BoneID b_id, NodeID n_id ){
            // ids exist
            assert( b_id >= 0 && b_id < bones.size() );
            assert( n_id >= 0 && n_id < nodes.size() );
            // bone is properly opened
            assert( nodes[bones[b_id].nodes.front()].type == SNT_Pole ||
                    nodes[bones[b_id].nodes.front()].type != SNT_Junction );
            // bone is not closed
            assert( nodes[bones[b_id].nodes.back()].type != SNT_Pole );
            assert( nodes[bones[b_id].nodes.back()].type != SNT_Junction );
            // node is not an endpoint
            assert( nodes[n_id].type != SNT_Pole );
            assert( nodes[n_id].type != SNT_Junction );
            // bone has at least 1 node
            assert( bones.back().nodes.size() > 0 );
            
            bones[b_id].nodes.push_back( n_id );
            buildNeighborhood( b_id );
            
        }
        
        void killCurrentBone( ){
            assert( bones.back().nodes.size() == 1 );
            bones.pop_back();
        }
        
        void shiftIDs( NodeID lastNodeID, BoneID lastBoneID ){
            for( int i = 0; i < bones.size(); ++i){
                bones[i].ID += lastBoneID;
                for( int ni = 0; i < bones[i].nodes.size(); ++ni ){
                    bones[i].nodes[ni] += lastNodeID;
                }
            }            
            for( int i = 0; i < nodes.size(); ++i){
                nodes[i].ID += lastNodeID;
                std::set<NodeID> new_neighbors;
                for( NodeID id : nodes[i].neighbors ){
                    new_neighbors.insert( id + lastNodeID );
                }
                nodes[i].neighbors = std::move( new_neighbors );
            }
        }
        
        bool sanityCheck(){
            
            bool ok = true;
            
            for( int i = 0; i < bones.size() && ok; ++i){
                ok = ok && ( nodes[bones[i].nodes.front()].type == SNT_Pole || nodes[bones[i].nodes.front()].type == SNT_Junction );
                ok = ok && ( nodes[bones[i].nodes.back()].type == SNT_Pole || nodes[bones[i].nodes.back()].type == SNT_Junction );
            }
            return ok;
        }
        
        public :
        std::vector< SkelNode > nodes;
        std::vector< SkelBone > bones;
        std::map< HMesh::VertexID, NodeID> poleToNode;
        std::map< HMesh::VertexID, NodeID> junctionSingularityToNode;
        
        void build( HMesh::Manifold& m, const std::set<HMesh::VertexID> &poleSet ){
            
            NodeID nodeID = 0;
            std::set< HMesh::VertexID> skelStarters;
            
            // build edge_info structure
            HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges(m);
            for(auto h: m.halfedges())
                if(edge_info[h].is_rib()) {
                    number_rib_edges(m, edge_info, h);
                    break;
                }
            
            // Find all Skel Leafs ( Poles )
            for( HMesh::VertexID v : poleSet ){
                skelStarters.insert( v );
                // save pole
                SkelNode n;
                n.ID        = nodeID;
                n.type      = SNT_Pole;
                n.pos       = m.pos( v );
                
                poleToNode[v] = n.ID;
                nodes.push_back( n );
                ++nodeID;
            }
            
            // Find all Skel Branching Nodes ( Complex Junction Loops )
            std::map< HMesh::HalfEdgeID, bool >     visited;
            std::map< HMesh::HalfEdgeID, NodeID >   junctionToNode;
            for( auto h : m.halfedges()){
                if( !edge_info[h].is_junction( )) { continue; } // skip all non-junction
                if ( visited[h] ){ continue; }                  // skip all visited
                
                junctionToNode[h]   = nodeID;
                double radius       = 0.0;
                CGLA::Vec3d centroid(0);
                std::set<HMesh::VertexID> vertices;
                std::queue< HMesh::HalfEdgeID > startHE;
                
                startHE.push( h );
                while( !startHE.empty() ){
                    HMesh::HalfEdgeID starter = startHE.front();
                    startHE.pop();
                    
                    HMesh::Walker w = m.walker( starter );
                    do{
                        assert(edge_info[w.halfedge()].is_junction());
                        assert(edge_info[w.opp().halfedge()].is_junction());
                        HMesh::HalfEdgeID curr      = w.halfedge(),
                        opp       = w.opp().halfedge();
                        HMesh::VertexID   curr_v    = w.vertex();
                        
                        visited[curr]        = true;
                        visited[opp]         = true;
                        junctionToNode[curr] = nodeID;
                        junctionToNode[opp]  = nodeID;
                        
                        if( vertices.count( curr_v )  == 0        ){
                            centroid += m.pos( curr_v );
                            vertices.insert( curr_v );
                        }
                        
                        if( valency( m, curr_v ) != 4 ){
                            junctionSingularityToNode[curr_v] = nodeID;
                            skelStarters.insert( curr_v );
                            HMesh::Walker sing = m.walker(curr_v);
                            for( ; !sing.full_circle(); sing = sing.circulate_vertex_ccw( )){
                                if( !visited[sing.halfedge()] && edge_info[sing.halfedge()].is_junction() ){
                                    startHE.push( sing.halfedge( ));
                                }
                            }
                        }
                        w = w.next().opp().next();
                    }while( w.halfedge() != starter );
                }
                centroid /= vertices.size();
                // calculate ball radius // use max radius
                for( HMesh::VertexID v : vertices ){
                    double length = ( centroid - m.pos(v) ).length();
                    if( length > radius ){ radius = length; }
                }
                assert( radius > 0.0 );
                // save pole
                SkelNode j;
                j.ID        = nodeID;
                j.type      = SNT_Junction;
                j.pos       = centroid;
                j.radius    = radius;
                nodes.push_back( j );
                ++nodeID;
            }
            
            // this should iterate through skelStarters and work differently if vertex is a pole or a junction
            //        for( HMesh::VertexID pole : poleSet ){
            std::map< HMesh::HalfEdgeID, NodeID >   ribToNode;
            for( HMesh::VertexID v : skelStarters ){
                
                std::queue< HMesh::HalfEdgeID > hes;
                NodeID startingNode;
                bool starter_is_pole = false;
                if( poleSet.count( v ) > 0 ){ // it is a pole
                    hes.push( m.walker(v).halfedge() );
                    assert( edge_info[hes.front()].is_spine( ));
                    startingNode    = poleToNode[v];
                    starter_is_pole = true;
                    
                }else{ // it is a singularity of the junction complex loop
                    assert( junctionSingularityToNode.count(v) > 0 );
                    startingNode = junctionSingularityToNode[v];
                    for( auto ws = m.walker( v ); !ws.full_circle(); ws = ws.circulate_vertex_cw() ){
                        if( edge_info[ws.halfedge()].is_spine( )){
                            hes.push( ws.halfedge( ));
                        }
                    }
                }
                
                while( !hes.empty()){
                    
                    auto he = hes.front();
                    hes.pop();
                    
                    HMesh::Walker w = m.walker( he );
                    openBone( startingNode );
                    HMesh::HalfEdgeID rib   = w.next().halfedge();
                    bool connectedToPole    = starter_is_pole;
                    
                    while( !edge_info[rib].is_junction() && ( poleSet.count( w.vertex()) == 0 )){
                        assert(edge_info[rib].is_rib());
                        assert(edge_info[w.halfedge()].is_spine());
                        
                        double radius       = 0.0;
                        CGLA::Vec3d centroid(0);
                        std::vector<HMesh::VertexID> loop_vertices;
                        HMesh::Walker loop_walker = m.walker( rib );
                        
                        if( ribToNode.count(rib) == 0){
                            for( ; !loop_walker.full_circle(); loop_walker = loop_walker.next().opp().next()){
                                
                                ribToNode[loop_walker.halfedge()]       = nodeID;
                                ribToNode[loop_walker.opp().halfedge()] = nodeID;
                                
                                centroid += m.pos( loop_walker.vertex( ));
                                loop_vertices.push_back( loop_walker.vertex( ));
                            }
                            centroid /= loop_vertices.size();
                            
                            for( HMesh::VertexID v : loop_vertices ){
                                double length = ( centroid - m.pos(v) ).length();
                                if( length > radius ){ radius = length; }
                            }
                            assert( radius > 0.0 );
                            SkelNode nl;
                            nl.ID       = nodeID;
                            nl.type     = SNT_Rib;
                            nl.pos      = centroid;
                            nl.radius   = connectedToPole ? ( centroid - nodes[startingNode].pos ).length() : radius;
                            nodes.push_back( nl );
                            addToCurrentBone( nl.ID );
                            ++nodeID;
                        }
                        else{ // the rib loop was already encountered. actually this situation should kill the current branch
                            // currently I'm deferring it to bone closure, but it should be killed here, performance gain!
                            killCurrentBone();
                            goto bone_killed;
                        }
                        
                        w   = w.next().opp().next();
                        rib = w.next().halfedge();
                        connectedToPole = false;
                    }
                    if( edge_info[rib].is_junction()){
                        assert(junctionToNode.count(rib) > 0 );
                        closeBone( junctionToNode[rib] );
                    }
                    else{
                        assert( poleSet.count(w.vertex()) > 0);
                        closeBone( poleToNode[w.vertex()] );
                    }
                    bone_killed : ;
                }
            }
        }
        
        void saveToFile( std::string path ){
            
            std::cout << " there are " << nodes.size() << " nodes  and " << bones.size() << " bones " << std::endl;
            
            std::ofstream t( path, std::ofstream::trunc );
            if( t.is_open()){
                t << "ID Cx Cy Cz RADIUS #NEIGHBORS NEIGHBORS_LIST" << std::endl << nodes.size() << std::endl;
                
                for( SkelNode& n : nodes ){
                    //                if( n.neighbors.size() == 0 ) { continue; }
                    t << n.ID << " " << n.pos[0] << " " << n.pos[1] << " " << n.pos[2] << " " << n.radius << " "
                    << n.neighbors.size();
                    for(NodeID neighbor : n.neighbors ){
                        t << " " << neighbor;
                    }
                    t << std::endl;
                    
                    // debug
                    std::string msg = ( n.type == SNT_Pole ? "Pole" : ( n.type == SNT_Junction ? "Junction" : "Rib") );
                    std::cout << n.ID << " ) is : " << msg << " # " << n.type <<  std::endl;
                }
                t.close();
            }
        }
        
        void reAlignIDs( const HMesh::VertexIDRemap& vremap ){

            std::map< HMesh::VertexID, NodeID> remapped;
            for( auto item : poleToNode ){
                assert(vremap.count( item.first ) > 0 );
                remapped[ vremap.at( item.first )] = item.second;
            }
            poleToNode = std::move( remapped );
        }
        
        // ther is copied by value
        void merge( const Skeleton &other, std::vector< std::pair< HMesh::VertexID, HMesh::VertexID > >& other_this_matches ){
            
            // this should connect the two skeletons at the match points
            // remove poles at each side
            // add a node in the middle with radius that is the mean between the radii of two opened tubes that will be glued.
            // extend the this-bone ended with one of match's pole in order to incorporate all the nodes of the other skeleton
            // connected to the other match's pole.
            // copy the other bones onto this skeleton
            
            // a lot of asserts
            
            // offset other
            // copy other's nodes and bones to this
            // fuse the bones
            std::vector< NodeID > this_pole_nodes, other_pole_nodes;
            for( auto item : other_this_matches ){
                
                assert( poleToNode.count( item.second ) > 0 );
                assert( other.poleToNode.count( item.first ) > 0 );
                
                this_pole_nodes.push_back( poleToNode[item.second] );
                other_pole_nodes.push_back( other.poleToNode.at( item.first ));
            }
            
            // add to this skeleton all other nodes except those involved in glueing
            // shifting their ID
            NodeID current_node_id = nodes.size() - 1;
            std::map< NodeID, NodeID > new_node_IDs;
            
            // copy the non-involved nodes
            for( const SkelNode& node : other.nodes ){
                if( std::find( other_pole_nodes.begin(), other_pole_nodes.end(), node.ID ) != other_pole_nodes.end()){ continue; }
                
                SkelNode new_node;
                new_node.ID     = ++current_node_id;
                new_node.radius = node.radius;
                new_node.pos    = node.pos;
                new_node.type   = node.type;
                
                nodes.push_back( new_node );
                new_node_IDs[node.ID] = new_node.ID;
            }
            // this make sense since those nodes will be merged
            for( int i = 0; i < this_pole_nodes.size(); ++i ){
                new_node_IDs[other_pole_nodes[i]] = this_pole_nodes[i];
            }
            
            std::vector< BoneID > other_bones_glued;
            
            // merge the glued poles and reopen their bones
            for( int i = 0; i< this_pole_nodes.size(); ++i ){
                const SkelNode& other_node = other.nodes[other_pole_nodes[i]];
                nodes[this_pole_nodes[i]].type = SNT_Rib;
                nodes[this_pole_nodes[i]].radius =
                    std::max( nodes[this_pole_nodes[i]].radius, other_node.radius );
                
                // copy the bone outgoing from other_pole_nodes[i]
                bool reverse            = other_node.ID == other.bones[other_node.boneID].nodes.back();
                size_t other_bone_size  = other.bones[other_node.boneID].nodes.size();
                BoneID curr_bone        = nodes[this_pole_nodes[i]].boneID;
                for( int bi = 1; bi < other.bones[other_node.boneID].nodes.size() - 1; ++i ){
                    int index = reverse ? bi : other_bone_size - bi - 1;

                    assert( new_node_IDs.count( other.bones[other_node.boneID].nodes[index] ) > 0 );
                    NodeID mapped_node = new_node_IDs[other.bones[other_node.boneID].nodes[index]];

                    assert( mapped_node >= 0 && mapped_node < nodes.size() );
                    addToBone( curr_bone, mapped_node );
                }
                if( reverse ){
                    closeBone( curr_bone, other.bones[other_node.boneID].nodes.front( ));
                }
                else{
                    closeBone( curr_bone, other.bones[other_node.boneID].nodes.back( ));
                }
                other_bones_glued.push_back( curr_bone );
            }
            
            for( const SkelBone& b : other.bones ){
                if( std::find( other_bones_glued.begin(), other_bones_glued.end(), b.ID ) != other_bones_glued.end()) { continue; }
                openBone(b.nodes.front());
                for( int bni = 0 ; bni < b.nodes.size() - 1; ++bni ){
                    assert( new_node_IDs.count(b.nodes[bni] > 0 ));
                    addToCurrentBone( new_node_IDs[b.nodes[bni]] );
                }
                closeBone( b.nodes.back( ));
            }
            assert( sanityCheck( ));
        }
    };

}

#endif
