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
        const NodeID start = bones.back().nodes.front();
        if( boneEndPoints[start].count(id) > 0 ){ // il bone esisteva già
            bones.pop_back();
        }else{
            if( boneEndPoints.count( start ) ==  0 )    { boneEndPoints[start] = std::set< NodeID >(); };
            if( boneEndPoints.count( id ) ==  0 )       { boneEndPoints[id]    = std::set< NodeID >(); };
            boneEndPoints[start].insert(id);
            boneEndPoints[id].insert( start );
        }

    }
    
    void addToCurrentBone( NodeID id ){
        assert( id >= 0 && id < nodes.size() );
        assert( nodes[id].type != SNT_Pole );
        assert( nodes[id].type != SNT_Junction );
        
        NodeID curr_last = bones.back().nodes.back();
        nodes[id].neighbors.insert( curr_last );
        nodes[curr_last].neighbors.insert( id );

        bones.back().nodes.push_back( id );
    }
    
    void removeCurrentBone(){
        assert(bones.back().nodes.size() > 0 );
        bones.pop_back();
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
    std::map< HMesh::VertexID,   NodeID> poleToNode;
    
    void build( HMesh::Manifold& m, PoleSet &poleSet ){
        
        std::map< HMesh::HalfEdgeID, NodeID> edgeToNode;
        
        NodeID nodeID = 0;
        
        // build edge_info structure
        HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges(m);
        for(auto h: m.halfedges())
            if(edge_info[h].is_rib()) {
                number_rib_edges(m, edge_info, h);
                break;
            }
        
        for( HMesh::VertexID v : poleSet ){
            // save pole
            SkelNode n;
            n.ID        = nodeID;
            n.type      = SNT_Pole;
            n.pos       = m.pos( v );
            
            poleToNode[v] = n.ID;
            nodes.push_back( n );
            ++nodeID;
        }
        
        for( HMesh::VertexID v : poleSet ){

            NodeID poleId = poleToNode[v];
            
            std::queue< HMesh::HalfEdgeID > starters;

            // save each outgoing edge from the pole
            for( HMesh::Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw()){
                starters.push( w.halfedge( ));
                assert( edge_info[w.halfedge()].is_spine( ));
            }
            assert( starters.size() == valency( m, v ));
            
            // for each outgoing pole of the edge
            while( !starters.empty( )){
                bool connectedToPole    = true;
                bool closeCurrentBone   = true;
                
                openBone( poleId );
                HMesh::Walker pole_to_pole = m.walker( starters.front() );
                // need to iterate until it arrives to a pole
                while( poleSet.count( pole_to_pole.vertex()) == 0){
                    
                    
                    assert(edge_info[pole_to_pole.halfedge()].is_spine());
                    
                    std::vector<CGLA::Vec3d>        loop_v_pos;
                    std::vector<HMesh::HalfEdgeID>  loop_he;
                    HMesh::HalfEdgeID rib = pole_to_pole.next().halfedge();
                    assert( edge_info[rib].is_rib() );
                    bool rib_is_junction = edge_info[rib].is_junction();
                    
                    // I can skip it because I've already calculated that loop' centroid + radius
                    if( edgeToNode.count( rib ) == 0 ){
                    
                        // calculate centroid and radius
                        HMesh::Walker w = m.walker( rib );
                        double radius   = 0.0;
                        CGLA::Vec3d centroid( 0 );
                        
                        while ( !w.full_circle( )){
                            // here nodeID is the one I will use for the node I'm creating
                            if( !rib_is_junction ){ // save correspondance halfedge->nodeId
                                edgeToNode[w.halfedge()]        = nodeID;
                                edgeToNode[w.opp().halfedge()]  = nodeID;
                            }
                            
                            CGLA::Vec3d pos = m.pos( w.vertex() );
                            loop_v_pos.push_back( pos );
                            centroid += pos;
                            w = w.next().opp().next();
                        }
                        centroid /= loop_v_pos.size();
                        
                        // calculate ball radius // use max radius
                        for( CGLA::Vec3d& v_pos : loop_v_pos ){
                            double length = ( centroid - v_pos ).length();
                            if( length > radius ){ radius = length; }
                        }
                        assert( radius > 0.0 );
                        
                        SkelNode nl;
                        nl.ID       = nodeID;
                        nl.type     = rib_is_junction ? SNT_Junction : SNT_Rib;
                        nl.pos      = centroid;
                        nl.radius   = connectedToPole ? ( centroid - nodes[poleId].pos ).length() : radius;
                        nodes.push_back( nl );
                        ++nodeID;
                        if( rib_is_junction ){ closeBone( nl.ID ); openBone( nl.ID ); }
                        else                 { addToCurrentBone( nl.ID ); }
                    }
                    // devo solo aggiungerlo al bone corrente
                    else{
                        addToCurrentBone( edgeToNode[rib] );
//                        removeCurrentBone();
//                        closeCurrentBone = false;
//                        break;
                    }
                    
                    connectedToPole = false;
                    pole_to_pole = pole_to_pole.next().opp().next();
                }
                
//                if( closeCurrentBone ){
                    HMesh::VertexID end_pole = pole_to_pole.vertex();
                    assert( poleSet.count( end_pole ) > 0 );
                    assert( is_pole( m, end_pole ));
                    assert( poleToNode.count(end_pole) > 0 );
                    
                    closeBone( poleToNode[end_pole] );
//                }
                
                starters.pop();
            }
        }
        assert(sanityCheck());
        
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
