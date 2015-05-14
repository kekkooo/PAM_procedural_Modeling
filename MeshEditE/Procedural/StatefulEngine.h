//
//  StatefulEngine.h
//  MeshEditE
//
//  Created by Francesco Usai on 17/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__StatefulEngine__
#define __MeshEditE__StatefulEngine__

#include <stdio.h>
#include <vector>
#include <map>
#include <random>
#include <set>

#include <GEL/HMesh/Manifold.h>

#include <GEL/Geometry/KDTree.h>

#include <GEL/CGLA/Vec3d.h>
#include <GEL/CGLA/Mat4x4d.h>

#include "MeshEditE/Procedural/Helpers/module_alignment.h"
#include "MeshEditE/Procedural/EngineHelpers/InfoContainers.h"
#include "MesheditE/Procedural/Module.h"
#include "MesheditE/Procedural/MainStructure.h"

namespace GEL_Geometry = Geometry;

namespace Procedural{
    namespace Engines{

    /************************************************
     * TYPE DECLARATIONS                            *
     ***********************************************/
typedef std::vector<HMesh::VertexID>                                    VertexList;
typedef std::map<HMesh::VertexID, CGLA::Vec3d>                          VertexPosMap;
typedef std::set<HMesh::VertexID>                                       VertexSet;
typedef std::map< HMesh::VertexID, HMesh::VertexID >                    VertexMatchMap;
typedef GEL_Geometry::KDTree< CGLA::Vec3d, HMesh::VertexID >            kD_Tree;
        
typedef std::pair< std::vector< Procedural::GraphMatch::Match>,
                                Procedural::GraphMatch::EdgeCost >      matchesAndCost;
typedef std::pair<HMesh::VertexID, double>                              IdDistPair;
typedef std::vector<IdDistPair>                                         IDsDistsVector;
typedef std::pair<double, GraphMatch::EdgeCost>                         ExtendedCost;
        
inline bool operator <( const ExtendedCost& l, const ExtendedCost& r)
{
    if( l.first < r.first && l.second < r.second )
    {
        return true;
    }
    else if( GraphMatch::is_equal( l.first, r.first ))
    {
        return ( l.second < r.second );
    }
    else if( l.second == r.second )
    {
        return ( l.first < r.first );
    }
    else { return false; }
}
        
// MATCH è ( polo_modulo, polo_host )

        
struct MatchInfoProxy{
    private:
        Procedural::Helpers::ModuleAlignment::match_info matchInfo;
        bool                                             isValid    = false;
    public:
        bool IsValid() { return isValid; }
        Procedural::Helpers::ModuleAlignment::match_info& getMatchInfo() {
            return  matchInfo;
        }
        void setMatchInfo( Procedural::Helpers::ModuleAlignment::match_info& mi )
        {
            matchInfo = std::move( mi );
            isValid = true;
        }
};
        
struct CandidateInfo{
    HMesh::VertexID id;
};
struct CandidateSet{
public:
    const VertexSet& getCandidates(){ return candidates; }
    const std::vector<HMesh::VertexID> getCandidateVector(){
        std::vector<HMesh::VertexID> v;
        for( VertexID c : candidates ){ v.push_back(c); }
        return v;
    }
    void insert( HMesh::VertexID id, CandidateInfo info ){
        candidates.insert( id );
        candidate_infos[id] = info;
    }
    CandidateInfo getInfo( HMesh::VertexID id ){
        assert( candidates.count(id) > 0 );
        return candidate_infos[id];
    }
    void clear(){
        candidates.clear();
        candidate_infos.clear();
    }
private:
    VertexSet candidates;
    HMesh::AttributeVector<CandidateInfo, HMesh::VertexID> candidate_infos;
};

/// Stateful Engine class
class StatefulEngine{
    
    enum DimensionalityConstraint { Constrained_1D, Constrained_2D, Constrained_3D };
    
    /************************************************
     * METHODS                                      *
     ***********************************************/
    public :
    static  StatefulEngine& getCurrentEngine();
            void            setHost( HMesh::Manifold &host );
            void            setModule( HMesh::Manifold &module );
            void            testMultipleTransformations( int no_tests, size_t no_glueings );
            void            glueModuleToHost();
            void            consolidate();
    
    
            void            applyRandomTransform();
            void            applyOptimalAlignment( );
            void            alignModuleNormalsToHost( );
            void            alignUsingBestMatch( );
            void            actualGlueing();
    
    private :
                            StatefulEngine();
                            StatefulEngine( StatefulEngine const& ) = delete;
            void operator   = (StatefulEngine const&)               = delete;

            void            buildRandomRotation( CGLA::Mat4x4d &t );
            // you can assume that the module's poles are vertices of its convex hull
            void            buildCollisionAvoidingTranslation(
                                const CGLA::Mat4x4d &rot, CGLA::Mat4x4d &tr );

            void            buildCollisionAvoidingRandomTransform( CGLA::Mat4x4d &t );

            /*  INHERITED FROM module_alignment */
            void            buildMainStructureKdTree();
            void            transformModulePoles( CGLA::Mat4x4d &t, VertexPosMap &new_pos );

            void            matchModuleToHost( Procedural::PoleInfoMap& poleInfoMap, VertexMatchMap& M_pole_to_H_vertex );
            bool            findSecondClosest( const HMesh::VertexID &pole, const PoleGeometryInfo &pgi,
                                               const HMesh::VertexID &closest, HMesh::VertexID &second_closest, VertexSet &assigned );
            void            addNecessaryPoles();
    
            void            buildTransformationList( std::vector< CGLA::Mat4x4d> &transformations );


    

    
    /************************************************
     * ATTRIBUTES                                   *
     ***********************************************/
private:

    HMesh::Manifold     *m;

    VertexSet           H_vertices,
                        M_vertices;

    
    Procedural::MainStructure*  mainStructure;
    Procedural::Module*         candidateModule;
    
    std::vector<Procedural::Module>
                        transformedModules;

    kD_Tree*            tree;
    bool                treeIsValid;
    
    Procedural::EngineHelpers::EdgeInfoContainer edge_info;
    
    MatchInfoProxy      best_match;
    
    std::mt19937_64     randomizer;
    double              last_x1, last_x2, last_x3;
    size_t              current_glueing_target;


};
        

    
}}
#endif /* defined(__MeshEditE__StatefulEngine__) */
