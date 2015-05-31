//
//  Module.h
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Module__
#define __MeshEditE__Module__

#include <stdio.h>

#include <map>
#include <vector>
#include <queue>
#include <set>
#include <fstream>
#include <iostream>

#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Mat4x4d.h>

#include <polarize.h>

namespace Procedural{
    
enum PoleLabel {
    NoneP    = 0,       // this is 00000000 it will match no value
    RedP     = 1,
    BlueP    = 2,
    GreenP   = 4,
    CyanP    = 8,
    MagentaP = 16,
    YellowP  = 32,
    PurpleP  = 64,
    AllP     = 127      // this is 11111111 it will match each value
};

typedef unsigned int Moduletype;
    
struct PoleGeometryInfo{
    unsigned int    valence;
    CGLA::Vec3d     pos;
    CGLA::Vec3d     normal;
};
    
struct PoleInfo{
    PoleGeometryInfo    geometry;
    Moduletype          moduleType  = 0;     // an identifier for the module's type
    PoleLabel           labels      = AllP;  // a label attached to a pole
    int                 age;                 // pole's starting age
    bool                isFree      = false;
};
    
typedef std::map<HMesh::VertexID, PoleInfo>             PoleInfoMap;
typedef std::vector< HMesh::VertexID>                   PoleList;
typedef std::set< HMesh::VertexID >                     PoleSet;

typedef size_t NodeID;
typedef size_t BoneID;

enum SkelNodeType{ SNT_Pole, SNT_Rib, SNT_Junction };
    
struct SkelNode{
    NodeID                  ID;
    SkelNodeType            type    = SNT_Rib;
    CGLA::Vec3d             pos;
    double                  radius  = 0.0;
    std::set<NodeID>        neighbors;
};

struct SkelBone{
    BoneID                  ID;
    std::vector<NodeID>   nodes;
};
    
struct Skeleton{
private:
    
    std::map< NodeID, std::set< NodeID > > boneEndPoints;
    
    void openBone( NodeID id ){
        assert( id >= 0 && id < nodes.size() );
        assert( nodes[id].type == SNT_Pole || nodes[id].type == SNT_Junction );
        SkelBone b;
        b.ID = bones.size();
        b.nodes.push_back( id );
        bones.push_back( b );
    }
    
    void closeBone( NodeID id ){
        assert( bones.back().nodes.size() > 0 );
        SkelNodeType t_front = nodes[bones.back().nodes.front()].type;
        assert( t_front == SNT_Pole || t_front == SNT_Junction );
        assert( id >= 0 && id < nodes.size() );
        assert( nodes[id].type == SNT_Pole || nodes[id].type == SNT_Junction );
        NodeID curr_last = bones.back().nodes.back();
        nodes[id].neighbors.insert( curr_last );
        nodes[curr_last].neighbors.insert( id );

        bones.back().nodes.push_back( id );
        const NodeID start = bones.back().nodes.front(),
                     end   = bones.back().nodes.back();
        
        if(( boneEndPoints[start].count(id) > 0 ) || ( boneEndPoints[end].count(id) > 0 ) ){ // bone already exists
            bones.pop_back();
        }else{
            if( boneEndPoints.count( start ) ==  0 )    { boneEndPoints[start] = std::set< NodeID >(); };
            if( boneEndPoints.count( id ) ==  0 )       { boneEndPoints[id]    = std::set< NodeID >(); };
            boneEndPoints[start].insert(id);
            boneEndPoints[id].insert( start );
            
            SkelBone& b = bones.back();

            for( int i = 1; i < b.nodes.size() - 1; ++i ){
                NodeID prev = b.nodes[i-1],
                       curr = b.nodes[i],
                       next = b.nodes[i+1];

                nodes[curr].neighbors.insert( prev );
                nodes[curr].neighbors.insert( next );
                nodes[prev].neighbors.insert( curr );
                nodes[next].neighbors.insert( curr );
            }
        }

    }
    
    void addToCurrentBone( NodeID id ){
        assert( id >= 0 && id < nodes.size() );
        assert( nodes[id].type != SNT_Pole );
        assert( nodes[id].type != SNT_Junction );
        
        bones.back().nodes.push_back( id );
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
    std::set< HMesh::VertexID> skelStarters;
    
    void build( HMesh::Manifold& m, PoleSet &poleSet ){
        
        NodeID nodeID = 0;
        
        // build edge_info structure
        HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges(m);
        for(auto h: m.halfedges())
            if(edge_info[h].is_rib()) {
                number_rib_edges(m, edge_info, h);
                break;
            }
        
        // Find all Skel Leafs ( Poles )
        for( HMesh::VertexID v : poleSet ){
            skelStarters.insert( v);
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
                    visited[w.halfedge()]               = true;
                    visited[w.opp().halfedge()]         = true;
                    junctionToNode[w.halfedge()]        = nodeID;
                    junctionToNode[w.opp().halfedge()]  = nodeID;
                    
                    if( vertices.count( w.vertex()) == 0        ){
                        centroid += m.pos( w.vertex( ));
                        vertices.insert( w.vertex( ));
                    }
                    
                    if( valency( m, w.vertex( )) != 4 ){
                        junctionSingularityToNode[w.vertex()] = nodeID;
                        skelStarters.insert( w.vertex() );
                        HMesh::Walker sing = m.walker(w.vertex());
                        for( ; !sing.full_circle(); sing = sing.circulate_vertex_ccw( )){
                            if( !visited[w.halfedge()] ){
                                startHE.push( w.halfedge( ));
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
        for( HMesh::VertexID v : skelStarters ){
            
            std::queue< HMesh::HalfEdgeID > hes;
            NodeID startingNode;
            bool starter_is_pole = false;
            if( poleSet.count(v) > 0 ){ // it is a pole
                hes.push( m.walker(v).halfedge() );
                assert( edge_info[hes.front()].is_spine( ));
                startingNode    = poleToNode[v];
                starter_is_pole = true;

            }else{
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
                    for( ; !loop_walker.full_circle(); loop_walker = loop_walker.next().opp().next()){
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

                    w   = w.next().opp().next();
                    rib = w.next().halfedge();
                    connectedToPole = false;
                }
                if( edge_info[rib].is_junction()){
                    closeBone( junctionToNode[rib] );
                }
                else{
                    assert( poleSet.count(w.vertex()) > 0);
                    closeBone( poleToNode[w.vertex()] );
                }
            }
        }
    }
    
    void saveToFile( std::string path ){
        
        std::cout << " there are " << nodes.size() << " nodes  and " << bones.size() << " bones " << std::endl;
        
        std::ofstream t( path, std::ofstream::trunc );
        if( t.is_open()){
            t << "ID Cx Cy Cz RADIUS #NEIGHBORS NEIGHBORS_LIST" << std::endl << nodes.size() << std::endl;
            
            for( SkelNode& n : nodes ){
                t << n.ID << " " << n.pos[0] << " " << n.pos[1] << " " << n.pos[2] << " " << n.radius << " "
                  << n.neighbors.size();
                for(NodeID neighbor : n.neighbors ){
                    t << " " << neighbor;
                }
                t << std::endl;
            }
            t.close();
        }
    }
    
    
    void glue( const Skeleton &other, const std::vector< std::pair< HMesh::VertexID, HMesh::VertexID > >& other_this_matches ){
        
        // this should connect the two skeletons at the match points
        // remove poles at each side
        // add a node in the middle with radius that is the mean between the radii of two opened tubes that will be glued.
        // extend the this-bone ended with one of match's pole in order to incorporate all the nodes of the other skeleton
        // connected to the other match's pole.
        // copy the other bones onto this skeleton
        
// TODOOOO
        // esportare gli scheletri nel cino-formato in modo da testarli.
        
    }
};
    

class Module{
    
public:
         // this will instantiate the internal manifold structure and pole info using obj_load
         Module () { m = NULL; }
         Module( std::string path, Moduletype mType );
         Module( HMesh::Manifold &manifold, Moduletype mType, size_t no_glueings = 1 );
    
        Module& getTransformedModule( const CGLA::Mat4x4d &T, bool transform_module = false );
        void reAlignIDs( HMesh::VertexIDRemap &remapper );
    
        const PoleInfo&    getPoleInfo( HMesh::VertexID p );
        bool isPole( HMesh::VertexID v );

        inline const PoleInfoMap& getPoleInfoMap(){ return poleInfoMap;}
private:
    void    BuildPoleInfo();            
    
/************************************************
* ATTRIBUTES                                   *
***********************************************/
public:
#warning those should be deleted from public
    HMesh::Manifold     *m;
    int                 no_of_glueings;
    PoleList            poleList;
    PoleSet             poleSet;

    CGLA::Vec3d         bsphere_center;
    double              bsphere_radius;
    
private :
    PoleInfoMap         poleInfoMap;
    Skeleton            *skeleton;

    };
}

#endif /* defined(__MeshEditE__Module__) */
