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

namespace GEL_Geometry = Geometry;

namespace Procedural{
    namespace Engines{

    /************************************************
     * TYPE DECLARATIONS                            *
     ***********************************************/
typedef std::vector<HMesh::VertexID>                                    VertexList;
typedef std::map<HMesh::VertexID, CGLA::Vec3d>                          VertexPosMap;
typedef std::set<HMesh::VertexID>                                       VertexSet;
typedef std::map< HMesh::VertexID, HMesh::VertexID >                    VertexMatch;
typedef GEL_Geometry::KDTree< CGLA::Vec3d, HMesh::VertexID >            kD_Tree;
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
                matchInfo = mi;
                isValid = true;
            }
        };

/// Stateful Engine class
class StatefulEngine{
    
    /************************************************
     * METHODS                                      *
     * devo aggiungere quelli relativi al branching *
     * e alla gestione del quanto di ogni cosa      *
     ***********************************************/
    public :
    static  StatefulEngine& getCurrentEngine();
            void            setHost( HMesh::Manifold &host );
            void            setModule( HMesh::Manifold &module );
            void            testMultipleTransformations( int no_tests, int no_glueings );
            void            glueModuleToHost();
            void            consolidate();
    
    private :
                            StatefulEngine();
                            StatefulEngine( StatefulEngine const& ) = delete;
            void operator   = (StatefulEngine const&)               = delete;

            void            buildRandomTransform( CGLA::Mat4x4d &t );
            // you can assume that the module's poles are vertices of its convex hull
            void            buildCollisionAvoidingTranslation(
                                const CGLA::Mat4x4d &rot, CGLA::Mat4x4d &tr );

            void            buildCollisionAvoidingRandomTransform( CGLA::Mat4x4d &t );

            /*  INHERITED FROM module_alignment */
            void            buildHostKdTree();
            void            transformModulePoles( CGLA::Mat4x4d &t, VertexPosMap &new_pos );

            /// this is to fix - look inside
            void            matchModuleToHost( VertexPosMap& module_poles_positions, VertexMatch& M_pole_to_H_vertex );
            void            findSecondClosest( const HMesh::VertexID &closest, HMesh::VertexID &second_closest,
                                               VertexSet &assigned);
            void            alignUsingBestMatch( );
            void            actualGlueing();


    

    
    /************************************************
     * ATTRIBUTES                                   *
     ***********************************************/
private:
    bool                treeIsValid;
    kD_Tree             tree;
    HMesh::Manifold     *m;
    
    VertexSet           H_vertices,
                        M_vertices;
    
    // VertexPosMap        M_original_pos_map;
    std::mt19937_64     randomizer;
//    Procedural::GraphMatch::match_info          best_match;
    MatchInfoProxy      best_match;


};
        

    
}}
#endif /* defined(__MeshEditE__StatefulEngine__) */
