//
//  HMeshParallelKit.cpp
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 07/01/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "HMeshParallelKit.h"

VertexIDBatches batch_vertices(const HMesh::Manifold& m) {
    VertexIDBatches vertex_ids(CORES);
    auto batch_size = m.no_vertices()/CORES;
    int cnt = 0;
    for_each_vertex(m, [&](HMesh::VertexID v) {
        vertex_ids[(cnt++/batch_size)%CORES].push_back(v);
    });
    return vertex_ids;
}
