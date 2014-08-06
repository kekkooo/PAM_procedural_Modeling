//
//  HMeshParallelKit.h
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 07/01/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__HMeshParallelKit__
#define __MeshEditE__HMeshParallelKit__

#include <iostream>
//
//  HMeshParallelKit.h
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 07/01/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef MeshEditE_HMeshParallelKit_h
#define MeshEditE_HMeshParallelKit_h

#include <thread>
#include <vector>
#include <GEL/HMesh/Manifold.h>

const int CORES = 16;


class range {
public:
    class iterator {
        friend class range;
    public:
        long int operator *() const { return i_; }
        const iterator &operator ++() { ++i_; return *this; }
        iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }
        
        bool operator ==(const iterator &other) const { return i_ == other.i_; }
        bool operator !=(const iterator &other) const { return i_ != other.i_; }
        
    protected:
        iterator(long int start) : i_ (start) { }
        
    private:
        unsigned long i_;
    };
    
    iterator begin() const { return begin_; }
    iterator end() const { return end_; }
    range(long int  begin, long int end) : begin_(begin), end_(end) {}
private:
    iterator begin_;
    iterator end_;
};


typedef std::vector<std::vector<HMesh::VertexID>> VertexIDBatches;

template<typename  T>
inline void for_each_vertex_parallel(int no_threads, const VertexIDBatches& batches, const T& f) {
    std::vector<std::thread> t_vec(no_threads);
    for(auto t : range(0, no_threads))
        t_vec[t] = std::thread(f, ref(batches[t]));
    for(auto t : range(0, no_threads))
        t_vec[t].join();
}

inline void for_each_vertex(HMesh::Manifold& m, std::function<void(HMesh::VertexID)> f) { for(auto v : m.vertices()) f(v); }
inline void for_each_vertex(const HMesh::Manifold& m, std::function<void(HMesh::VertexID)> f) { for(auto v : m.vertices()) f(v); }

VertexIDBatches batch_vertices(const HMesh::Manifold& m);

#endif


#endif /* defined(__MeshEditE__HMeshParallelKit__) */
