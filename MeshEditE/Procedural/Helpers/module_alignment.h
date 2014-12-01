//
//  module_alignment.h
//  MeshEditE
//
//  Created by Francesco Usai on 02/11/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__module_alignment__
#define __MeshEditE__module_alignment__

#include <stdio.h>

#include <map>
#include <set>

#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/CGLA/Mat4x4d.h>

#include <MeshEditE/Procedural/Helpers/svd_alignment.h>
#include <MeshEditE/Procedural/Matches/graph_match.h>

namespace GEL_Geometry = Geometry;

namespace Procedural {
    namespace Helpers{
        namespace ModuleAlignment{

#define TEST_PATH "/Users/francescousai/Documents/Dottorato/Visiting/Results/add/"
// REMEMBER - Procedural::GraphMatch::Match is a pair ( module's pole, host's vertex )

            
            
typedef GEL_Geometry::KDTree< CGLA::Vec3d, HMesh::VertexID >            kd_tree;
typedef std::map< HMesh::VertexID, HMesh::VertexID >                    vertex_match;
typedef std::pair< std::vector< Procedural::GraphMatch::Match>,
                                Procedural::GraphMatch::EdgeCost >      matches_and_cost;
typedef std::pair<HMesh::VertexID, double>                              ID_and_dist;
typedef std::vector<ID_and_dist>                                        IDs_and_dists;
struct match_info{
    Procedural::GraphMatch::EdgeCost            cost;
    std::vector<Procedural::GraphMatch::Match>  matches;
    CGLA::Mat4x4d                               random_transform;
};
        
void AddModule(                 HMesh::Manifold &host, HMesh::Manifold &module, size_t no_glueings,
                                std::vector<Procedural::GraphMatch::Match> &matches );
void build_manifold_kdtree(     HMesh::Manifold& m, std::set< HMesh::VertexID > &selected, kd_tree &tree);
CGLA::Mat4x4d transform_module_poles(
                                HMesh::Manifold &host, const std::set<HMesh::VertexID> &host_vs,
                                const std::set<HMesh::VertexID> &module_vs,
                                std::map< HMesh::VertexID, CGLA::Vec3d > &new_pos );
void generate_random_transform( HMesh::Manifold &host, HMesh::Manifold &module, CGLA::Mat4x4d &t );
/// this is to fix - look inside
void match_module_to_host(      HMesh::Manifold &m, kd_tree &tree, std::map< HMesh::VertexID, CGLA::Vec3d > &module_poles_positions,
                                vertex_match &pole_to_host_vertex );
            
Procedural::GraphMatch::EdgeCost
     get_glueings(              HMesh::Manifold &host, HMesh::Manifold &module, vertex_match &pole_to_host_vertex,
                                size_t no_glueings, std::vector< HMesh::VertexID > &host_p,
                                std::vector< HMesh::VertexID > &module_p );
void apply_optimal_alignment(   HMesh::Manifold &m, const std::set<HMesh::VertexID> &module_vs, match_info &choosen_match_info );
void best_configuration(        std::vector<matches_and_cost> matches_and_costs, matches_and_cost &choosen_cost );
void add_necessary_poles(       HMesh::Manifold &m, std::vector<Procedural::GraphMatch::Match> &matches,
                                std::set<HMesh::VertexID> &host_vs, std::set<HMesh::VertexID> &module_vs );
void find_second_closest(       const HMesh::Manifold &m, const kd_tree &tree, const HMesh::VertexID &closest,
                                HMesh::VertexID &second_closest, std::set<HMesh::VertexID> &assigned);
void align_module_normals_to_host(
                                HMesh::Manifold &m, std::set<HMesh::VertexID> &module_IDs, std::vector<Procedural::GraphMatch::Match> &matches );
void glue_matches(              HMesh::Manifold &m, std::vector<Procedural::GraphMatch::Match> &matches );

}}}

#endif /* defined(__MeshEditE__module_alignment__) */
