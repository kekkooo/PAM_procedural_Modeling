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
#include <GEL/HMesh/Manifold.h>
#include <GEL/Geometry/KDTree.h>
#include <GEL/CGLA/Mat4x4d.h>
#include <map>

namespace Procedural {
    namespace Helpers{
        namespace ModuleAlignment{

typedef Geometry::KDTree< CGLA::Vec3d, HMesh::VertexID >    kd_tree;
typedef std::map< HMesh::VertexID, HMesh::VertexID >        vertex_match;
        
void AddModule(                 HMesh::Manifold &host,  HMesh::Manifold &module, size_t no_glueings );
void build_manifold_kdtree(     HMesh::Manifold& m,    std::vector< HMesh::VertexID > &selected, kd_tree &tree);
void transform_module_poles(    HMesh::Manifold &host, HMesh::Manifold &module,
                                std::map< HMesh::VertexID, CGLA::Vec3d > &new_pos );
void generate_random_transform( HMesh::Manifold &host, HMesh::Manifold &module, CGLA::Mat4x4d &t );
/// this is to fix - look inside
void match_module_to_host(      HMesh::Manifold &host, HMesh::Manifold &module, kd_tree &tree,
                                vertex_match &pole_to_host_vertex );
void get_glueings(              HMesh::Manifold &host, HMesh::Manifold &module, vertex_match &pole_to_host_vertex,
                                int no_glueings, std::vector< HMesh::VertexID > selected_poles );
        
        
        
}}}

#endif /* defined(__MeshEditE__module_alignment__) */
