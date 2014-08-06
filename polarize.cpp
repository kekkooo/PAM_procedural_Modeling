//
//  polarize.cpp
//  GEL
//
//  Created by J. Andreas Bærentzen on 18/03/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//
#include <set>
#include <queue>
#include <iomanip>
#include <thread>

#include <unistd.h>

#include <GEL/HMesh/harmonics.h>
#include "polarize.h"
#include <GEL/CGLA/statistics.h>
#include <GEL/CGLA/eigensolution.h>
#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/Mat4x4d.h>
#include <GEL/CGLA/Vec3i.h>
#include <GEL/LinAlg/LapackFunc.h>
#include <GEL/HMesh/triangulate.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/HMesh/curvature.h>
#include <GEL/HMesh/quadric_simplify.h>
#include <GEL/HMesh/mesh_optimization.h>
#include <GEL/HMesh/smooth.h>
#include <GEL/HMesh/dual.h>

#include <GEL/Geometry/KDTree.h>
#include <GEL/Geometry/QEM.h>
#include <GEL/Geometry/GridAlgorithm.h>
#include <GEL/Geometry/build_bbtree.h>
#include <GEL/Geometry/RGrid.h>
#include <GEL/Geometry/TrilinFilter.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>

#include <GEL/LinAlg/LapackFunc.h>

#include "LogMap.h"
#include "HMeshParallelKit.h"
#include <GEL/GLGraphics/MeshEditor.h>
#include "polarize.h"

using namespace CGLA;
using namespace Geometry;
using namespace std;
using namespace HMesh;
using namespace LinAlg;
using namespace GLGraphics;

bool is_pole(Manifold& m, VertexID v)
{
    Walker w = m.walker(v);
    for(; !w.full_circle(); w = w.circulate_vertex_ccw())
        if(w.face() == InvalidFaceID || no_edges(m, w.face()) != 3)
            return false;
    return true;
}


const Vec3f& get_color(int i)
{
    gel_srand(0);
    static Vec3f ctable[100000];
    static bool was_here;
    if(!was_here)
    {
        
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        gel_rand();
        was_here = true;
        ctable[0] = Vec3f(0);
        for(int j=1;j<100000;++j) {
            ctable[j] = Vec3f(gel_rand(),gel_rand(),gel_rand());
            ctable[j] /= GEL_RAND_MAX;
        }
    }
    return ctable[i%100000];
}

inline bool same_level(double a, double b) {return abs(a-b) < 1e-10;}

enum VertexStatus {
    RAW, COOKED
};

typedef KDTree<Vec3d, int> PtTree;

PtTree extract_poles(Manifold& m, VertexAttributeVector<double>& fun)
{
    PtTree tree;
    for(auto vid: m.vertices())
    {
        double v0 = fun[vid];
        bool maximum = true;
        bool minimum = true;
        for(Walker w = m.walker(vid); !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            if(fun[w.vertex()]>v0)
                maximum = false;
            if(fun[w.vertex()]<v0)
                minimum = false;
        }
        if(minimum || maximum) {
            tree.insert(m.pos(vid), 0);
            if(minimum)
                cout << setprecision(20) << "Minimum at: " << fun[vid] << endl;
            if(maximum)
                cout << "Maximum at: " << fun[vid] << endl;
        }
    }
    tree.build();
    return tree;
}

void compute_edge_weights(const Manifold& m, HalfEdgeAttributeVector<double>& edge_weights, FaceAttributeVector<int>& included)
{
    edge_weights = HalfEdgeAttributeVector<double>(m.allocated_halfedges(), 0);
    for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
        if(included[*f])
        {
            for(Walker wv = m.walker(*f); !wv.full_circle(); wv = wv.circulate_face_ccw())
            {
                HalfEdgeID h = wv.halfedge();
                Vec3d p1(m.pos(wv.vertex()));
                Vec3d p2(m.pos(wv.next().vertex()));
                Vec3d p0(m.pos(wv.opp().vertex()));
                double ang = acos(min(1.0, max(-1.0, dot(normalize(p1-p0), normalize(p2-p0)))));
                double ang_opp = acos(min(1.0, max(-1.0, dot(normalize(p2-p1), normalize(p0-p1)))));
                double l = (p1-p0).length();
                //                double a = acos(min(1.0, max(-1.0, dot(normalize(p0-p2), normalize(p1-p2)))));
                //                double w = max(0.0000001,0.5/tan(a));
                //                edge_weights[h]  += w;
                //                edge_weights[wv.opp().halfedge()]  += w;
                edge_weights[h] += tan(ang/2) / l;
                edge_weights[wv.opp().halfedge()] += tan(ang_opp/2) / l;
            }
        }
}


template<typename T>
void smooth_fun(const Manifold& m,
                const VertexAttributeVector<int>& nailed,
                VertexAttributeVector<T>& fun, int iter)
{
    auto batches = batch_vertices(m);
    HalfEdgeAttributeVector<double> edge_weights;
    FaceAttributeVector<int> included(m.allocated_faces(),1);
    compute_edge_weights(m,edge_weights, included);
    bool ignore_nailing = false;
    auto new_fun = fun;
    auto f = [&](const vector<VertexID>& vids) {
        for(VertexID v: vids)
            if(!nailed[v] || ignore_nailing)
            {
                double w_sum = 0;
                new_fun[v] = T(0);
                circulate_vertex_ccw(m, v, [&](Walker wv) {
                    double w = edge_weights[wv.halfedge()];
                    new_fun[v] += w * fun[wv.vertex()];
                    w_sum += w;
                });
                new_fun[v] /= w_sum;
                new_fun[v] = 0.5 * new_fun[v] + 0.5 * fun[v];
            }
    };
    for(int i = 0; i < iter-5; ++i)
    {
        for_each_vertex_parallel(CORES, batches, f);
        swap(fun,new_fun);
    }
    ignore_nailing = true;
    
    for(int i = 0; i < 5; ++i) {
        for_each_vertex_parallel(CORES, batches, f);
        swap(fun,new_fun);
    }
}

template void smooth_fun<double>(const Manifold& m,
                                 const VertexAttributeVector<int>& nailed,
                                 VertexAttributeVector<double>& fun, int iter);
template void smooth_fun<Vec3d>(const Manifold& m,
                                 const VertexAttributeVector<int>& nailed,
                                 VertexAttributeVector<Vec3d>& fun, int iter);

void shortest_edge_triangulate_face(Manifold& m, FaceID f0, const VertexAttributeVector<int>& ls_id)
{
    queue<FaceID> face_queue;
    
    face_queue.push(f0);
    
    while(!face_queue.empty())
    {
        FaceID f = face_queue.front();
        face_queue.pop();
        
        // Create a vector of vertices.
        vector<VertexID> verts;
        for(Walker w = m.walker(f); !w.full_circle(); w = w.circulate_face_ccw())
        {
            FaceID fa = w.face();
            FaceID fb = f;
            assert(fa==fb);
            verts.push_back(w.vertex());
        }
        // If there are just three we are done.
        if(verts.size() == 3) continue;
        
        // Find vertex pairs that may be connected.
        vector<pair<int,int> > vpairs;
        const int N = verts.size();
        for(int i = 0; i < N - 2; ++i){
            for(int j = i + 2; j < N; ++j){
                if(verts[i] != verts[j] &&
                   !connected(m, verts[i], verts[j]) &&
                   (ls_id[verts[i]] != ls_id[verts[j]] || ls_id[verts[i]]==0))
                    vpairs.push_back(pair<int,int>(i, j));
            }
        }
        if(vpairs.empty()){
            cout << "Warning: could not triangulate a " << verts.size() << "  face." << endl;
            
            for(int i=0;i<verts.size();++i)
                cout << verts[i] << " " << ls_id[verts[i]] << ",";
            cout << endl;
                
            
            continue;
        }
        
        /* For all vertex pairs, find the edge lengths. Combine the
         vertices forming the shortest edge. */
        
        float min_len=FLT_MAX;
        int min_k = -1;
        for(size_t k = 0; k < vpairs.size(); ++k){
            int i = vpairs[k].first;
            int j = vpairs[k].second;
            float len = valency(m,verts[i]) + valency(m,verts[j]);//length(m.pos(verts[i]) - m.pos(verts[j]));
            //            cout << m.pos(verts[i]) << m.pos(verts[j]) << len << endl;
            if(len<min_len){
                min_len = len;
                min_k = k;
            }
        }
        assert(min_k != -1);
        
        if(min_k>=0)
        {
            // Split faces along edge whose midpoint is closest to isovalue
            int i = vpairs[min_k].first;
            int j = vpairs[min_k].second;
            //           cout << i << " " << j << " " << min_k << " " << vpairs.size() << endl;
            FaceID f_new = m.split_face_by_edge(f, verts[i], verts[j]);
            
            if(no_edges(m, f)>3)
                face_queue.push(f);
            if(no_edges(m, f_new)>3)
                face_queue.push(f_new);
        }
        else
            cout << "failed to triangle. Suspect NaN vertex positions!!" << endl;
        
    }
}

class EdgeLengthEnergy: public EnergyFun
{
public:
    double delta_energy(const Manifold& m, HalfEdgeID h) const
    {
        Walker w = m.walker(h);
		
        VertexID hv = w.vertex();
		VertexID hov = w.opp().vertex();
		VertexID hnv = w.next().vertex();
		VertexID honv = w.opp().next().vertex();
		
		Vec3d va(m.pos(hv));
		Vec3d vb(m.pos(hov));
		Vec3d vc(m.pos(hnv));
		Vec3d vd(m.pos(honv));
        int val_ab = valency(m, hv) + valency(m, hov);
        int val_cd = 2+valency(m, hnv) + valency(m, honv);
        double ratio = 1;// val_ab/(val_ab+val_cd);
        //        return -((1-ratio)*sqr_length(vc-vd)-ratio*sqr_length(va-vb));
//        return (val_cd/val_ab)*length(vc-vd)-(val_ab/val_cd)*length(va-vb);
                return length(vc-vd)-length(va-vb);
    }
};

void flip_edges(Manifold& m, int id, const VertexAttributeVector<int>& ls_id, const VertexAttributeVector<double>& fun)
{
    vector<HalfEdgeID> hvec;
    for(HalfEdgeID hid : m.halfedges())
    {
        Walker w = m.walker(hid);
        if(ls_id[w.vertex()] == id &&
           ls_id[w.vertex()] != ls_id[w.opp().vertex()] &&
           no_edges(m, w.face())==3 && no_edges(m, w.opp().face())==3)
            hvec.push_back(hid);
    }
    
    EdgeLengthEnergy ele;
    ValencyEnergy vae;
    MinAngleEnergy mae(-1);
    
    bool did_work = false;
    int iter=0;
    double avglen = 0;
    for(auto h: m.halfedges())
        avglen += length(m, h);
    avglen /= m.no_halfedges();
    
    double T = 100 * avglen;
    do {
        T *= 0.9;
        did_work = false;
        random_shuffle(begin(hvec), end(hvec));
        for(auto hid : hvec)
        {
            Walker w = m.walker(hid);
            double delta_e = ele.delta_energy(m, hid);
            double rand_num = gel_rand()/double(GEL_RAND_MAX);
            if((ls_id[w.next().vertex()] != ls_id[w.opp().next().vertex()]) &&
               (delta_e<0 || rand_num < exp(-delta_e/T))) {
                m.flip_edge(hid);
                did_work = true;
            }
        }
    }
    while (did_work && ++iter<1000);
    cout << "edges " << hvec.size() << " id " << id;
    if(hvec.size())
        cout << " fun " << fun[m.walker(hvec[0]).vertex()];
    cout << " flip iter" << iter << endl;
}


void connect_touched(Manifold& m, const VertexAttributeVector<int>& touched)
{
    vector<tuple<FaceID,VertexID,VertexID>> connect_vec;
    for(FaceID fid : m.faces()) {
        VertexID v0 = InvalidVertexID;
        VertexID v1 = InvalidVertexID;
        for(Walker w = m.walker(fid);!w.full_circle(); w = w.next())
        {
            if(touched[w.vertex()]) {
//                assert(v1 == InvalidVertexID);
                if(v0 == InvalidVertexID)
                    v0 = w.vertex();
                else
                    v1 = w.vertex();
            }
        }
        if(v1 != InvalidVertexID)
            connect_vec.push_back(make_tuple(fid, v0, v1));
    }
    
     for(auto t: connect_vec)
        m.split_face_by_edge(get<0>(t), get<1>(t), get<2>(t));

}

vector<VertexID> cut_mesh(Manifold& m, VertexAttributeVector<double>& fun, double cut_val,
                          VertexAttributeVector<VertexStatus>& status)
{
    cout << "cutting @ " << cut_val << endl;
    
    vector<VertexID> new_verts;
    vector<HalfEdgeID> hidvec;
    VertexAttributeVector<int> touched(m.allocated_vertices(), 0);
    for(HalfEdgeID hid : m.halfedges())
    {
        Walker w = m.walker(hid);
        double bval = fun[w.vertex()];
        double aval = fun[w.opp().vertex()];
        if(aval<bval && aval < cut_val && cut_val <= bval)
        {
            double t = (cut_val-aval)/(bval-aval);
            Vec3d p = t*m.pos(w.vertex()) + (1.0-t)*m.pos(w.opp().vertex());
            VertexID vnew = m.split_edge(hid);
            m.pos(vnew) = p;
            status[vnew] = COOKED;
            fun[vnew] = cut_val;
            touched[vnew] = 1;
            new_verts.push_back(vnew);
        }
    }
    
    connect_touched(m, touched);
    
    return new_verts;
}


struct SkeletalLaplacian
{
    Manifold& m;
    typedef vector<pair<VertexID, double>> NeighborHood;
    VertexAttributeVector<NeighborHood> nbrs;
    VertexAttributeVector<Vec3d> tangents;
    VertexAttributeVector<double> fun;
    VertexAttributeVector<Mat3x3d> proj_tensor;
  
    double delta_f = 0.25;
    SkeletalLaplacian(Manifold& _m, const VertexAttributeVector<double>& _fun): m(_m), fun(_fun) {
        for(VertexID v0: m.vertices())
        {
            double f0 = fun[v0];
            double f_min = f0 - delta_f;
            double f_max = f0 + delta_f;
            FaceAttributeVector<int> faces_visited(m.allocated_faces(), 0);
            VertexAttributeVector<int> vertices_visisted(m.allocated_vertices(), 0);
            queue<FaceID> fq;
            
            circulate_vertex_ccw(m, v0, [&](FaceID f) {fq.push(f);});
            
            while(!fq.empty())
            {
                FaceID f = fq.front();
                fq.pop();
                
                if(!faces_visited[f])
                {
                    faces_visited[f] = 1;
                    circulate_face_ccw(m, f, [&](Walker& w){
                        VertexID v = w.vertex();
                        if(!vertices_visisted[v] &&
                           f_min <= fun[v] && fun[v] < f_max) {
                            nbrs[v0].push_back(make_pair(v, fun[v] - fun[v0]));
                            vertices_visisted[v] = 1;
                        }
                        double aval = fun[w.vertex()];
                        double bval = fun[w.opp().vertex()];
                        if(! ((aval < f_min && bval < f_min) || (aval > f_max && bval > f_max))) {
                            FaceID fo = w.opp().face();
                            if(!faces_visited[fo])
                                fq.push(fo);
                        }
                    });
                }
            }
            
            size_t N = nbrs[v0].size();
            cout << "nbrs :" << N << endl;
            
            CMatrix U(N,4);
            CVector b(N);
            CVector x(N);

            Vec3d mean(0);
            for(int i=0;i<N;++i)
            {
                VertexID v = nbrs[v0][i].first;
                Vec3d p = m.pos(v);
                mean += p;
                U[i][0] = p[0];
                U[i][1] = p[1];
                U[i][2] = p[2];
                U[i][3] = 1;
                b[i] = fun[v];
            }
            CMatrix UT = U.Transposed();
            CMatrix D = Inverted(UT * U)*UT;
            x = D * b;
            tangents[v0] = normalize(Vec3d(x[0], x[1], x[2]));
            mean /= N;
            
            vector<Vec3d> pts;
            for(int i=0;i<N;++i)
            {
                VertexID v = nbrs[v0][i].first;
                Vec3d p = m.pos(v)-mean;
                Vec3d t = tangents[v0];
                p -= t * dot(p,t);
                pts.push_back(p);
            }
            
            
            Mat3x3d C,Q,L;
            covariance(pts, C);
            int n = power_eigensolution(C, Q, L);
            double w = (L[0][0]/L[1][1]);
            cout << "n = " << n << " w = " << w << endl;
            outer_product(w*Q[1], w*Q[1], proj_tensor[v0]);
            
            double wsum = 0;
            bool bigger= false;
            bool smaller=false;
            for(int i=0;i<N;++i)
            {
                if(nbrs[v0][i].second >0) bigger=true;
                if(nbrs[v0][i].second <0) smaller=true;
                
                nbrs[v0][i].second = exp(-abs(nbrs[v0][i].second));
                wsum += nbrs[v0][i].second;
            }
            for(int i=0;i<N;++i)
                nbrs[v0][i].second /= wsum;
            
            if(!(bigger && smaller)) nbrs[v0].resize(0);
        }
    }
    
    void smooth_proj_dirs() {
        for(VertexID v0: m.vertices())
        {
            Mat3x3d mat = 30*proj_tensor[v0];
            double w_sum = 30;
            for(auto n: nbrs[v0]) {
                VertexID v = n.first;
                double w = exp(-sqr(fun[v0]-fun[v]));
                mat += w*proj_tensor[v];
                w_sum += w;
            }
            proj_tensor[v0] = mat/w_sum;
        }
    
    }
    
    Vec3d get_proj_dir(VertexID v0) {
        Mat3x3d Q,L;
        power_eigensolution(proj_tensor[v0], Q, L);
        return Q[0];
    }
    
    
    VertexAttributeVector<Vec3d> apply(const VertexAttributeVector<Vec3d>& vec_i) {
        VertexAttributeVector<Vec3d> vec_o;
        for(VertexID v: m.vertices()) {
            Vec3d L = Vec3d(0);
            Vec3d p0 = vec_i[v];
            for(auto& n: nbrs[v]) {
                L += n.second * (vec_i[n.first] - p0);
            }
            vec_o[v] = L; // nbrs[v].size();
        }
        return  vec_o;
    }
    
};


void skeletonize(Manifold& m, const VertexAttributeVector<double>& fun, double w)
{
    SkeletalLaplacian L(m, fun);
    
    VertexAttributeVector<Vec3d> l = L.apply(m.positions_attribute_vector());

    for(VertexID v: m.vertices())
        m.pos(v) += w*l[v];
}

void skeleton_line_field(Manifold& m, const VertexAttributeVector<double>& fun, VertexAttributeVector<Vec3d>& lin)
{
    SkeletalLaplacian L(m, fun);
    for(int i=0;i<100;++i)
        L.smooth_proj_dirs();
    for(VertexID v: m.vertices())
        lin[v] = L.get_proj_dir(v);
}




void label_connected_components(Manifold&m, int& cur_id, const vector<VertexID>& vertices,
                                const VertexAttributeVector<double>& fun,
                                const VertexAttributeVector<VertexStatus>& status,
                                VertexAttributeVector<int>& ls_id)
{
    for(auto vid : vertices)
    {
        if(ls_id[vid]>0)
            continue;
        
        queue<VertexID> vq;
        vq.push(vid);
        while(!vq.empty())
        {
            VertexID v = vq.front();
            vq.pop();
            ls_id[v] = cur_id;
            for(Walker w = m.walker(v); !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                if(status[w.vertex()]==COOKED && same_level(fun[v], fun[w.vertex()])) {
                    if(ls_id[w.vertex()]==0)
                        vq.push(w.vertex());
                    else if(ls_id[w.vertex()] != cur_id)
                        cout << "previously labeled vertex encountered" << endl;
                }
            }
        }
        ++cur_id;
    }
}

void separate_edges(Manifold& m, int id,
                    VertexAttributeVector<double>& fun,
                    VertexAttributeVector<int>& ls_id)
{
    vector<HalfEdgeID> edges_to_split;

    HalfEdgeID h0 = InvalidHalfEdgeID;
    for(auto h: m.halfedges())
    {
        Walker w = m.walker(h);
        if(ls_id[w.vertex()]==id &&
           ls_id[w.opp().vertex()] == id &&
           abs(fun[w.vertex()]) > abs(fun[w.next().vertex()]))
        {
            h0 = h;
            break;
        }
    }
    if(h0 == InvalidHalfEdgeID)
    {
        cout << " Could not separate edges; id = " << id << endl;
        return;
    }
    
    
    Walker w = m.walker(h0);
    
    while(!w.full_circle()) {
        if(ls_id[w.vertex()] == id)
            w = w.next();
        else
        {
            edges_to_split.push_back(w.opp().halfedge());
            w = w.opp().next();
        }
    }
    
    VertexAttributeVector<int> touched(m.allocated_vertices(),0);
    vector<VertexID> kill_vertices;
    for (auto h : edges_to_split)
    {
        VertexID v = m.walker(h).vertex();
        VertexID vnew = m.split_edge(h);
        fun[vnew] = fun[v];
        m.pos(vnew) = m.pos(v);
        ls_id[vnew] = id;
        touched[vnew] = 1;

        ls_id[v] = 0;
        kill_vertices.push_back(v);
    }


    connect_touched(m, touched);

    for(auto vid: kill_vertices)
        {

            FaceID fid = m.merge_one_ring(vid);

            if(fid != InvalidFaceID)
                shortest_edge_triangulate_face(m, fid, ls_id);
        }
}

typedef vector<pair<double, VertexID>> LSCurve;

LSCurve trace_ls_curve(Manifold& m, const VertexAttributeVector<int>& ls_id, int id)
{
    LSCurve lsc;
    for(auto vid : m.vertices()) {
        if(m.in_use(vid) && ls_id[vid]==id) {
            lsc.push_back(make_pair(0, vid));
            VertexID v0 = vid;
            Walker w = m.walker(vid);
            do {
                for(; ls_id[w.vertex()] != id ; w=w.circulate_vertex_ccw());
                double t = lsc.back().first;
                t += length(m,w.halfedge());
                lsc.push_back(make_pair(t, w.vertex()));
                w = w.opp().prev().opp();
            }
            while (lsc.back().second != v0);
            break;
        }
    }
    return lsc;
}

void smooth_loop_vertices(Manifold& m, const VertexAttributeVector<int>& ls_id, int id, int max_iter=12)
{
    auto new_pos = m.positions_attribute_vector();
    LSCurve lsc = trace_ls_curve(m, ls_id, id);
    LSCurve lsc0 = lsc;
    
    if(lsc0.size()==0)
        return;
    
    for(int iter=0; iter < max_iter; ++ iter)
    {
        LSCurve lsc_new = lsc;
        for(int i=1;i<lsc.size()-1;++i)
        {
            lsc_new[i].first = 0.5*(lsc[i].first + 0.5*(lsc[i+1].first+lsc[i-1].first));
        }
        lsc = lsc_new;
    }
    for(int i=1;i<lsc.size()-1;++i)
    {
        double t = lsc[i].first;
        for(int j=0;j<lsc0.size();++j)
        {
            double t0 = lsc0[j].first;
            double t1 = lsc0[j+1].first;
            if(t0<t && t<t1)
            {
                double weight = (t-t0)/(t1-t0);
                new_pos[lsc[i].second] =
                (1-weight)*m.pos(lsc0[j].second) +
                weight * m.pos(lsc0[j+1].second);
                break;
            }
        }
        
    }
    
    m.positions_attribute_vector() = new_pos;
}

void remove_v5_vertices(Manifold& m, const VertexAttributeVector<int>& ls_id)
{
    queue<HalfEdgeID> Q;
    for(auto h: m.halfedges())
    {
        Walker w = m.walker(h);
        
        int id = ls_id[w.vertex()];
        if(id==ls_id[w.opp().vertex()])
        {
                vector<HalfEdgeID> hit_list;
                w=w.circulate_vertex_ccw();
                for(; ls_id[w.vertex()] != id; w=w.circulate_vertex_ccw())
                    hit_list.push_back(w.halfedge());
                
                if(hit_list.size()>1)
                {
                    sort(begin(hit_list),end(hit_list),
                         [&m](HalfEdgeID a, HalfEdgeID b){return length(m,a)<length(m, b);});
                    for(int i=1;i<hit_list.size();++i)
                    {
                        HalfEdgeID hid = hit_list[i];
                        Walker wn = m.walker(hid);
                        if(valency(m,wn.vertex())==4)
                            Q.push(wn.next().opp().next().halfedge());
                        m.merge_faces(wn.face(), hid);
                    }
                }
        }
    }
    
    while (!Q.empty()) {
        HalfEdgeID h = Q.front();
        Q.pop();
        
        if(m.in_use(h)) {
            Walker w= m.walker(h);
            
            if(valency(m, w.vertex()) != 3)
                Q.push(w.next().opp().next().halfedge());
            m.merge_faces(w.face(), h);
        }
    }
}

void remove_wedges(Manifold& m, int id, const VertexAttributeVector<int>& ls_id, const VertexAttributeVector<double>& fun)
{
    vector<HalfEdgeID> hidvec;
    for(auto hid : m.halfedges())
    {
        Walker w = m.walker(hid);
        if(ls_id[w.vertex()] == id &&
           ls_id[w.vertex()] == ls_id[w.opp().vertex()] &&
           abs(fun[w.next().vertex()])<=abs(fun[w.vertex()]) &&
           no_edges(m, w.face())==3 &&
           precond_collapse_edge(m, hid))
            hidvec.push_back(hid);
    }
    
    random_shuffle(begin(hidvec), end(hidvec));
    
    for(auto hid: hidvec)
        if(m.in_use(hid) &&
           precond_collapse_edge(m, hid))
            m.collapse_edge(hid, true);
    
    
}

void remove_vertices(Manifold& m, const VertexAttributeVector<VertexStatus>& status,
                     const VertexAttributeVector<int>& ls_id)
{
    for(auto f: m.faces())
        if(no_edges(m, f)>4)
            shortest_edge_triangulate_face(m, f, ls_id);
    
    vector<VertexID> vertices_to_push;
    for(auto vid: m.vertices())
        if(status[vid]==RAW)
            vertices_to_push.push_back(vid);
    
    random_shuffle(begin(vertices_to_push), end(vertices_to_push));
    
    for(auto vid : vertices_to_push)
    {
        FaceID fid = m.merge_one_ring(vid);
        if(fid != InvalidFaceID)
            shortest_edge_triangulate_face(m, fid, ls_id);
    }
}

//void polarize_mesh(Manifold& m_in, VertexAttributeVector<double>& fun, double vmin, double vmax, const int divisions)
//{
//
//    auto f = [](Manifold& m, VertexAttributeVector<double>& fun, double vmin, double vmax, int divisions) {
//        parallel_work.lock(); // ---------------------
//        PtTree pole_tree = extract_poles(m, fun);
//        
//        VertexAttributeVector<VertexStatus> status(m.allocated_vertices(), RAW);
//        VertexAttributeVector<int> ls_id(m.allocated_vertices(), 0);
//        int cur_id = 1;
//        
//        auto f = [&](double cut_val) {
//            vector<VertexID> cut_verts = cut_mesh(m, fun,cut_val, status);
//            label_connected_components(m, cur_id, cut_verts, fun, status, ls_id);
//        };
//        
//        f(0);
//        parallel_work.unlock(); // ++++++++++++++++++++
////        std::this_thread::sleep_for(std::chrono::milliseconds(10));
//        
//        int start_id = cur_id;
//        double interval = (vmax-vmin)/(divisions+1);
//        for(int i=1;i<divisions;++i) {
//            parallel_work.lock(); // ---------------------
//            f(i*interval);
//            f(-i*interval);
//            parallel_work.unlock(); // ++++++++++++++++++++
//            std::this_thread::sleep_for(std::chrono::milliseconds(5));
//
//        }
//        
//        parallel_work.lock(); // ---------------------
//        remove_vertices(m, status, ls_id);
//        parallel_work.unlock(); // ++++++++++++++++++++
//        std::this_thread::sleep_for(std::chrono::milliseconds(5));
//
//
//        for(int id=1;id<cur_id;++id) {
//            parallel_work.lock(); // ---------------------
//            smooth_loop_vertices(m, ls_id, id,10);
//            flip_edges(m, id, ls_id, fun);
//            parallel_work.unlock(); // ++++++++++++++++++++
//            std::this_thread::sleep_for(std::chrono::milliseconds(5));
//        }
//        for(int id=start_id;id<cur_id;++id)
//        {
//            parallel_work.lock(); // ---------------------
//            separate_edges(m, id, fun, ls_id);
//            remove_wedges(m, id, ls_id, fun);
//            parallel_work.unlock(); // ++++++++++++++++++++
//            std::this_thread::sleep_for(std::chrono::milliseconds(5));
//
//
//            parallel_work.lock(); // ---------------------
//            smooth_loop_vertices(m, ls_id, id,30);
//            parallel_work.unlock(); // ++++++++++++++++++++
//            std::this_thread::sleep_for(std::chrono::milliseconds(5));
//
//        }
//        parallel_work.lock(); // ---------------------
//        remove_v5_vertices(m, ls_id);
//
//        parallel_work.unlock(); // ++++++++++++++++++++
////        std::this_thread::sleep_for(std::chrono::milliseconds(10));
//
//        parallel_work.lock(); // ---------------------
//        dual(m);
//
//        parallel_work.unlock(); // ++++++++++++++++++++
////        std::this_thread::sleep_for(std::chrono::milliseconds(10));
//
//        parallel_work.lock(); // ---------------------
//        for(auto vid : m.vertices())
//            if(is_pole(m, vid))
//            {
//                Vec3d k;
//                int v;
//                double dist=DBL_MAX;
//                pole_tree.closest_point(m.pos(vid), dist, k, v);
//                m.pos(vid) = k;
//            }
//        m.cleanup();
//        parallel_work.unlock(); // ++++++++++++++++++++
////        std::this_thread::sleep_for(std::chrono::milliseconds(10));
//
//    };
//    
//    static thread t;
//    if(t.joinable())
//        t.join();
//    t = thread(f,ref(m_in), ref(fun), vmin, vmax, divisions);
//    t.join();
//}


void polarize_mesh(Manifold& m, VertexAttributeVector<double>& fun, double vmin, double vmax, const int divisions)
{
    PtTree pole_tree = extract_poles(m, fun);
    
    VertexAttributeVector<VertexStatus> status(m.allocated_vertices(), RAW);
    VertexAttributeVector<int> ls_id(m.allocated_vertices(), 0);
    int cur_id = 1;
    
    auto f = [&](double cut_val) {
        vector<VertexID> cut_verts = cut_mesh(m, fun,cut_val, status);
        label_connected_components(m, cur_id, cut_verts, fun, status, ls_id);
    };
    
    if(vmin<0)
        f(0);
    
    int start_id = cur_id;
    double interval = (vmax-vmin)/divisions;
    for(int i=1;i<divisions+1;++i) {
        f(i*interval);
        
        if(vmin<0)
            f(-i*interval);
    }
    
    remove_vertices(m, status, ls_id);
    for(int id=1;id<cur_id;++id) {
        smooth_loop_vertices(m, ls_id, id,10);
        flip_edges(m, id, ls_id, fun);
    }

    for(int id=start_id;id<cur_id;++id)
    {
        separate_edges(m, id, fun, ls_id);
        remove_wedges(m, id, ls_id, fun);
        smooth_loop_vertices(m, ls_id, id,30);
    }
    remove_v5_vertices(m, ls_id);
    dual(m);
    for(auto vid : m.vertices())
        if(is_pole(m, vid))
        {
            Vec3d k;
            int v;
            double dist=DBL_MAX;
            pole_tree.closest_point(m.pos(vid), dist, k, v);
            m.pos(vid) = k;
        }
    m.cleanup();
}

Region trace_region(Manifold& m, HalfEdgeID hid)
{
    Walker w = m.walker(hid);
    if(no_edges(m, w.face()) == 4 ||
       no_edges(m, w.opp().face()) == 4)
    {
        double a_max=0;
        double a_sum=0;
        vector<HalfEdgeID> hvec;
        bool abort=false;
        for(int i=0;i<2;++i)
        {
            Walker wlk = w;
            if(i==0) {
                hvec.push_back(wlk.halfedge());
                double l =length(m,wlk.halfedge());
                a_sum += l;
                a_max = max(a_max, l);
            }
            else
                wlk = w.opp();
            
            while(no_edges(m, wlk.face())==4)
            {
                wlk = wlk.next().next().opp();
                if(wlk.halfedge() == w.halfedge()) {
                    i = 2;
                    break;
                }
                hvec.push_back(wlk.halfedge());
                double l =length(m,wlk.halfedge());
                a_sum += l;
                a_max = max(a_max, l);
            }
            if (no_edges(m, wlk.face())==3 && valency(m,wlk.next().vertex())<=3)
                abort = true;
        }
        return make_pair(abort?-1:/*(a_sum/hvec.size())*/a_max,hvec);
    }
    return make_pair(-1.0, vector<HalfEdgeID>(0));
}

void collapse_region(Manifold& m, Region& region)
{
    for(HalfEdgeID h: region.second)
    {
        if(m.in_use(h) && precond_collapse_edge(m, h))
            m.collapse_edge(h,true);
        else cout << "Collapse failed ..." << endl;
    }
}

void refine_region(Manifold& m, Region& region)
{
    VertexAttributeVector<int> touched(m.allocated_faces(),0);
    for(HalfEdgeID h: region.second)
    {
        Walker w = m.walker(h);
        if(no_edges(m, w.face())==3)
            touched[w.next().vertex()] = 1;
        if(no_edges(m,w.opp().face())==3)
            touched[w.opp().next().vertex()] = 1;
        if(m.in_use(h)) {
            VertexID v = m.split_edge(h);
            touched[v] = 1;
        }
    }
    connect_touched(m, touched);
}

void simplify_polar_mesh(Manifold& m, double frac, int max_iter)
{
    for(int iter=0;iter<max_iter;++iter) {
        int N_orig = m.no_vertices();
        HalfEdgeAttributeVector<int> touched(m.allocated_halfedges(),0);
        
        vector<Region> work_vector;
        for(auto hid : m.halfedges())
            if(!touched[hid])
            {
                auto region = trace_region(m, hid);
                if (region.first>0)
                    work_vector.push_back(region);
                for(auto h: region.second)
                {
                    Walker w = m.walker(h);
                    touched[w.halfedge()]=1;
                    touched[w.opp().halfedge()]=1;
                }
            }
        
        auto itr = min_element(begin(work_vector), end(work_vector));
        if(itr != end(work_vector))
            collapse_region(m,*itr);
    }
    m.cleanup();
    
}


HalfEdgeAttributeVector<EdgeInfo> label_PAM_edges(Manifold& m)
{
    HalfEdgeAttributeVector<EdgeInfo> edge_info(m.allocated_halfedges());
    queue<HalfEdgeID> hq;
    for(VertexID vid : m.vertices())
        if(is_pole(m, vid))
            circulate_vertex_ccw(m, vid, [&](Walker w) {
                edge_info[w.halfedge()] = EdgeInfo(SPINE, 0);
                edge_info[w.opp().halfedge()] = EdgeInfo(SPINE, 0);
                hq.push(w.opp().halfedge());
            });
    
    while(!hq.empty())
    {
        HalfEdgeID h = hq.front();
        Walker w = m.walker(h);
        hq.pop();
        bool is_spine = edge_info[h].edge_type == SPINE;
        for(;!w.full_circle(); w=w.circulate_vertex_ccw(),is_spine = !is_spine)
            if(edge_info[w.halfedge()].edge_type == UNKNOWN)
            {
                EdgeInfo ei = is_spine ? EdgeInfo(SPINE,0) : EdgeInfo(RIB,0);
                edge_info[w.halfedge()] = ei;
                edge_info[w.opp().halfedge()] = ei;
                hq.push(w.opp().halfedge());
            }
        
    }
    return edge_info;
}

int number_rib_edges(Manifold& m,  HalfEdgeAttributeVector<EdgeInfo>& edge_info, HalfEdgeID h)
{
    int rib_id=0;
    
    queue<HalfEdgeID> Q;
    Q.push(h);
    
    while(!Q.empty())
    {
        HalfEdgeID hid = Q.front();
        Q.pop();
        if(edge_info[hid].edge_type == RIB && edge_info[hid].id==0)
        {
            ++rib_id;
            bool branch = false;
            vector<Walker> w_vec;
            queue<Walker> wq;
            Walker w = m.walker(hid);
            wq.push(w);
            wq.push(w.opp());
            while(!wq.empty())
            {
                Walker w=wq.front();
                wq.pop();
                edge_info[w.halfedge()] = EdgeInfo(RIB, rib_id);
                edge_info[w.opp().halfedge()] = EdgeInfo(RIB, rib_id);
                if(valency(m, w.vertex())>4)
                    branch = true;
                assert(valency(m, w.vertex())%2 == 0);
                w_vec.push_back(w);
                for( w = w.next().opp().next(); edge_info[w.halfedge()].id == 0; w = w.opp().next().opp().next()) {
                    wq.push(w);
                    Q.push(w.next().next().opp().halfedge());
                    Q.push(w.next().opp().next().opp().halfedge());
                }
                
                for(auto w : w_vec)
                    if (branch)
                    {
                        edge_info[w.halfedge()].edge_type = RIB_JUNCTION;
                        edge_info[w.opp().halfedge()].edge_type = RIB_JUNCTION;
                    }
                
            }
        }
    }
    return rib_id+1;
}

FaceAttributeVector<int> segment_faces(Manifold& m, const HalfEdgeAttributeVector<EdgeInfo>& edge_info)
{
    FaceAttributeVector<int> face_segment(m.allocated_faces(), 0);
    queue<FaceID> fq;
    int seg_no=0;
    vector<vector<FaceID>> segments(1);
    for(auto fid: m.faces())
    {
        if(face_segment[fid]==0)
        {
            ++seg_no;
            segments.push_back(vector<FaceID>());
            fq.push(fid);
            face_segment[fid] = seg_no;
            segments[seg_no].push_back(fid);
            while(!fq.empty())
            {
                FaceID fid = fq.front();
                fq.pop();
                circulate_face_ccw(m, fid, [&](Walker w){
                    if(edge_info[w.halfedge()].edge_type != RIB_JUNCTION && face_segment[w.opp().face()]==0) {
                        fq.push(w.opp().face());
                        face_segment[w.opp().face()] = seg_no;
                        segments[seg_no].push_back(w.opp().face());
                    }
                });
            }
        }
    }
    for(int i=1;i<=seg_no;++i) {
        cout << "select -r ";
        for(auto f : segments[i])
            cout << ".f[" << f << "] ";
        cout << ";" << endl;
        cout << "sets -e -forceElement lambert" << (i+1) << "SG;" << endl;
    }
    return face_segment;
}


void polar_segment(Manifold& m, bool show_segments)
{
    m.cleanup();
    
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges(m);
    for(auto h: m.halfedges())
        if(edge_info[h].is_rib()) {
            number_rib_edges(m, edge_info, h);
            break;
        }
    FaceAttributeVector<int> face_segment;
    if (show_segments) {
        face_segment = segment_faces(m, edge_info);
    }
    for(HalfEdgeID hid : m.halfedges()) {
        if (edge_info[hid].edge_type == SPINE) {
            Vec3f col;
            if (show_segments)
                col = 0.9*get_color(face_segment[m.walker(hid).face()]);
            else col = Vec3f(0.8,0,0);
            if(DebugRenderer::edge_colors[hid] != Vec3f(0,1,0))
                DebugRenderer::edge_colors[hid] = col;
        }
        else {
            const Vec3f col_normal(0,0,0.6);
            const Vec3f col_branch(.25,.25,1);
            Vec3f col = edge_info[hid].edge_type == RIB_JUNCTION ? col_branch : col_normal;
            Walker w = m.walker(hid);
            if (show_segments && edge_info[hid].edge_type != RIB_JUNCTION)
                col = 0.9 * get_color(face_segment[w.face()]);
            if(DebugRenderer::edge_colors[w.halfedge()] != Vec3f(0,1,0)) {
                DebugRenderer::edge_colors[w.halfedge()] = col;
                DebugRenderer::edge_colors[w.opp().halfedge()] = col;
            }
            DebugRenderer::vertex_colors[w.vertex()] = col;
            DebugRenderer::vertex_colors[w.opp().vertex()] = col;
        }
        
    }
    for(FaceID fid : m.faces()) {
        Vec3f col(1);
        if (show_segments) col = get_color(face_segment[fid]);
        DebugRenderer::face_colors[fid] = col;
    }
    for(VertexID vid : m.vertices())
        if(is_pole(m, vid))
            DebugRenderer::vertex_colors[vid] = Vec3f(0.6,0,0);
}


struct SkeletonNode {
    vector<int> children = vector<int>(0);
    int parent = 0;
    bool in_use = true;
    int num = 0;
    Vec3d pos = Vec3d(0);
};

VertexAttributeVector<map<int,float>> skinning;

void show_skin(Manifold& m, int j)
{
    for(VertexID v: m.vertices())
    {
        if(skinning[v][j]>0)
            DebugRenderer::vertex_colors[v] = Vec3f(1,0,0) * skinning[v][j];
        else
            DebugRenderer::vertex_colors[v] = Vec3f(0,0,0.3);
    }
}

void Skinning(Manifold& m, int no_ribs,
              const HalfEdgeAttributeVector<EdgeInfo>& edge_info,
              const vector<SkeletonNode>& skel_nodes)
{
//    // Linear numbering of all vertices
//    int k = 0;
//    VertexAttributeVector<int> vertex_numbers(m.allocated_vertices(),0);
//    for(VertexID v: m.vertices())
//        vertex_numbers[v] = k++;
    
    // Skinning of all vertices except poles
    
    VertexAttributeVector<int> rib_id;
    for(HalfEdgeID h: m.halfedges())
    {
        if((edge_info[h].edge_type == RIB ||
            edge_info[h].edge_type == RIB_JUNCTION))
        {
            auto w = m.walker(h);
            VertexID v = w.vertex();
            int id = edge_info[h].id;
            rib_id[v] = id;
            
            // Skinning a vertex on a rib loop used as joint
            if(skel_nodes[id].in_use)
                skinning[v][id] = 1.0;
            // Skinning vertex on other rib loop
            else {
                vector<pair<float, int>> w_n; // weight node combos.
                queue<pair<int,Walker>> q; // Spine queue
                
                // Push outgoing spine edges onto queue
                int vval = valency(m, v)/2;
                w = w.next();
                for(int i=0;i<vval;++i) {
                    q.push(make_pair(1, w));
                    w = w.opp().next().opp().next();
                }
                
                // Traverse along all spine paths (branch where needed) until a rib loop
                // With associated joint is encountered.
                while(!q.empty()) {
                    int steps = q.front().first;
                    Walker wlk = q.front().second;
                    q.pop();
                    int id = edge_info[wlk.next().halfedge()].id;
                    if(skel_nodes[id].in_use)
                        w_n.push_back(make_pair(1.0/(steps+1), id));
                    else {
                        VertexID vid = wlk.vertex();
                        wlk = wlk.next().opp().next();
                        int vval = valency(m, vid)/2-1;
                         for(int i=0;i<vval;++i) {
                            q.push(make_pair(steps+1,wlk));
                            wlk=wlk.opp().next().opp().next();
                        }
                    }
                }
                double w_sum=0;
                for(auto& x: w_n)
                    w_sum += x.first;
                for(auto& x: w_n)
                    skinning[v][x.second] = x.first/w_sum;
            }
            
        }

    }
    
    // Skinning pole
    for(VertexID v: m.vertices())
    {
        if(is_pole(m, v)) {
            Walker w = m.walker(v);
            int id =edge_info[w.next().halfedge()].id;
            skinning[v][id] = 1.0;
        }
    }
    
    
    // Smoothing the skinning
    for(int iter=0;iter<20;++iter)
        for(VertexID v: m.vertices())
        {
            circulate_vertex_ccw(m, v, [&](VertexID vn){
                for(const auto& x: skinning[vn])
                    skinning[v][x.first] += x.second;
            });
            double wsum = 0;
            for(const auto& x: skinning[v])
                wsum += x.second;
            for(auto& x: skinning[v])
                x.second /= wsum;
        }
   
    // Map to Maya numbers
    map<int,int> jnm;
    {int i=0;
        jnm[1]=i++; jnm[5]=i++; jnm[45]=i++; jnm[106]=i++; jnm[111]=i++; jnm[133]=i++; jnm[175]=i++; jnm[176]=i++; jnm[182]=i++; jnm[6]=i++; jnm[10]=i++; jnm[50]=i++; jnm[56]=i++; jnm[7]=i++; jnm[39]=i++; jnm[51]=i++; jnm[92]=i++; jnm[102]=i++; jnm[107]=i++; jnm[122]=i++; jnm[156]=i++; jnm[157]=i++; jnm[165]=i++; jnm[32]=i++; jnm[53]=i++; jnm[69]=i++; jnm[104]=i++; jnm[131]=i++; jnm[167]=i++; jnm[178]=i++; jnm[168]=i++; jnm[169]=i++; jnm[179]=i++; jnm[170]=i++; jnm[60]=i++; jnm[80]=i++; jnm[110]=i++; jnm[125]=i++; jnm[163]=i++; jnm[171]=i++; jnm[180]=i++; jnm[172]=i++; jnm[181]=i++; jnm[173]=i++; jnm[68]=i++; jnm[98]=i++; jnm[103]=i++; jnm[142]=i++; jnm[177]=i++; jnm[159]=i++; jnm[160]=i;




    }
    
    int k = 0;
    for(VertexID v: m.vertices()) {
        vector<float> weights(jnm.size(),0);
        for(auto x: skinning[v])
            weights[jnm[x.first]] = x.second;
        cout << "setAttr -s " << (jnm.size()) << " \".wl[" << k << "].w[0:" << jnm.size()-1 << "]\"";
        for(auto x: weights) {
            cout << " " << x;
        }
        cout << ";" <<endl;

        k++;
    }
}

void skeleton_retract(Manifold& m, double frac, VertexID v)
{
    HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges(m);
    int no_ribs;
    
    if (v == InvalidVertexID) {
        for(auto h: m.halfedges())
            if(edge_info[h].is_rib()) {
                no_ribs = number_rib_edges(m, edge_info, h);
                break;
            }
    }
    else {
        HalfEdgeID h;
        circulate_vertex_ccw(m, v, [&](HalfEdgeID hv){if(edge_info[hv].is_rib()) h = hv;});
        no_ribs = number_rib_edges(m, edge_info, h);
    }

    vector<SkeletonNode> skel_nodes(no_ribs);
    skel_nodes[1].parent = -1;
    for(HalfEdgeID hid : m.halfedges())
        if(edge_info[hid].edge_type == RIB ||
           edge_info[hid].edge_type == RIB_JUNCTION)
        {
            int id = edge_info[hid].id;
            SkeletonNode& node = skel_nodes[id];
            Walker w = m.walker(hid);
            node.pos += m.pos(w.vertex());
            node.num += 1;
            if(node.parent==0) {
                bool r1 = edge_info[w.next().next().halfedge()].is_rib();
                bool r2 = edge_info[w.opp().next().next().halfedge()].is_rib();
                if(r1 && r2)
                    node.parent = min(edge_info[w.next().next().halfedge()].id,
                                     edge_info[w.opp().next().next().halfedge()].id);
                else if(r1)
                    node.parent= edge_info[w.next().next().halfedge()].id;
                else
                    node.parent = edge_info[w.opp().next().next().halfedge()].id;
                
                skel_nodes[node.parent].children.push_back(id);
            }
        }
    for(int i=0;i<no_ribs; ++i)
        skel_nodes[i].pos /= skel_nodes[i].num;
    for(HalfEdgeID hid : m.halfedges())
        if(edge_info[hid].is_rib())
        {
            m.pos(m.walker(hid).vertex()) *= (1-frac);
            m.pos(m.walker(hid).vertex()) += (frac)*skel_nodes[edge_info[hid].id].pos;

        }
    
    auto remove_skeleton_node = [&](int n) {
        SkeletonNode& node = skel_nodes[n];
        if(node.in_use && node.children.size()>0)
        {
            int p = node.parent;
            for(int j=0; j<node.children.size();++j)
                skel_nodes[node.children[j]].parent = p;
            SkeletonNode& parent = skel_nodes[p];
            parent.children.erase(find(begin(parent.children), end(parent.children), n));
            parent.children.insert(end(parent.children), begin(node.children), end(node.children));
            node.in_use = false;
        }
    };
    
    for(int iter=0;iter<4;++iter)
        for(int i=2;i<no_ribs;++i)
            if(skel_nodes[i].children.size()>1)
                remove_skeleton_node(i);
    
    for(int iter=0;iter<14;++iter)
        for(int i=2;i<no_ribs;++i) {
            SkeletonNode& node = skel_nodes[i];
            if(node.children.size()==1)
            {
                Vec3d v0 = node.pos-skel_nodes[node.children[0]].pos;
                Vec3d v1 = skel_nodes[node.parent].pos-node.pos;
                v0.normalize();
                v1.normalize();
                if(dot(v0,v1)>0.96)
                    remove_skeleton_node(i);
            }
        }
    
    for(HalfEdgeID h: m.halfedges()) {
        if((edge_info[h].edge_type == RIB ||
           edge_info[h].edge_type == RIB_JUNCTION))
        {
            auto w = m.walker(h);
            VertexID v = w.vertex();
            int id = edge_info[h].id;
            if(id == 1)
                DebugRenderer::vertex_colors[v] = Vec3f(1,0.7,.7);
            else if(skel_nodes[id].in_use)
                DebugRenderer::vertex_colors[v] = Vec3f(1);
            else
                DebugRenderer::vertex_colors[v] = Vec3f(0);
            DebugRenderer::face_colors[w.face()] = Vec3f(0.3);
        }
        DebugRenderer::edge_colors[h] = Vec3f(0.1);
    }

    for(int i=1;i<no_ribs;++i) {
        SkeletonNode& node = skel_nodes[i];
        if(node.in_use) {
            Vec3d p = node.pos;
            if(node.parent > 0)
                cout << "select -r joint" << node.parent << ";" << endl;
            cout << "joint -n joint" << i << " "
            << "-p " << p[0] << " " << p[1] << " " << p[2] << ";" << endl;
            if(node.parent > 0)
                cout << "joint -e -zso -oj xyz -sao yup joint" << node.parent << ";" << endl;
        }
    }
    
   Skinning(m, no_ribs, edge_info, skel_nodes);
    
}


bool pole_quadrangulation(Manifold& m)
{
    vector<VertexID> poles;
    for(VertexID v: m.vertices())
        if(is_pole(m,v))
            poles.push_back(v);
    
    vector<VertexID> new_vertices;
    
    auto fan_compress = [&](VertexID pole, vector<HalfEdgeID>& hvec)
    {
        while(hvec.size()>=2)
        {
            for(auto h : hvec)
                m.split_edge(h);
            
            vector<HalfEdgeID> clpse_hvec;
            for(int i=0;i<hvec.size()-1;++i)
            {
                Walker w = m.walker(hvec[i]);
                FaceID f = w.face();
                VertexID v0 = w.opp().vertex();
                VertexID v1 = w.next().next().vertex();
                FaceID fn = m.split_face_by_edge(f, v0, v1);
                clpse_hvec.push_back(m.walker(fn).halfedge());
            }
            
            VertexID v_new;
            for(auto h: clpse_hvec)
            {
                v_new = m.walker(h).vertex();
                if(precond_collapse_edge(m, h))
                    m.collapse_edge(h);
            }
            new_vertices.push_back(v_new);
            m.pos(v_new) = m.pos(pole);
            hvec = vector<HalfEdgeID>(hvec.begin()+1, hvec.end()-1);
        }
        
        if(hvec.size()==1)
            m.merge_faces(m.walker(hvec[0]).face(), hvec[0]);
        
    };
    
    for(VertexID pole : poles) {
        int cw, ccw, n = valency(m, pole) - 2;
        if(n % 2 == 0)
            cw = ccw = n/2;
        else {
            ccw = n/2+1;
            cw = n/2;
        }
        
        HalfEdgeID h0=m.walker(pole).halfedge();
        double l0=length(m,h0);
        circulate_vertex_ccw(m, pole, [&](HalfEdgeID h){
            double l = length(m,h);
            if(l<l0) {
                h0 = h;
                l0 = l;
            }
        });
        
        vector<HalfEdgeID> ccw_vec;
        Walker w = m.walker(h0).circulate_vertex_ccw();
        for(int i=0;i<ccw;++i) {
            ccw_vec.push_back(w.halfedge());
            w = w.circulate_vertex_ccw();
        }
        fan_compress(pole, ccw_vec);
        vector<HalfEdgeID> cw_vec;
        
        w = m.walker(h0).circulate_vertex_cw();
        for (int i=0; i<cw; ++i) {
            cw_vec.push_back(w.halfedge());
            w = w.circulate_vertex_cw();
        }
        reverse(begin(cw_vec), end(cw_vec));
        fan_compress(pole, cw_vec);
    }
    
    double avg_hlen = 0;
    for(auto h : m.halfedges())
        avg_hlen += length(m, h);
    avg_hlen /= m.no_halfedges();
    
    
    VertexAttributeVector<Vec3d> new_pos = m.positions_attribute_vector();
    for(int iter=0;iter<50;++iter) {
        VertexAttributeVector<Vec3d> lap;
        for(VertexID vn: m.vertices()) {
            Vec3d avg(0);
            double wsum = 0;
            circulate_vertex_ccw(m, vn, [&](VertexID vnn){
                double wgt = exp(-sqr_length(m.pos(vn)-m.pos(vnn))/sqr(12*avg_hlen));
                avg+=wgt*m.pos(vnn);
                wsum += wgt;
            });
            lap[vn] = avg/wsum - m.pos(vn);
        }
       
        for(VertexID vn: new_vertices) {
            Vec3d sqlap(0);
            double w=0;
            double wsum = 0;
            circulate_vertex_ccw(m, vn, [&](VertexID vnn) {
                double wgt = exp(-sqr_length(m.pos(vn)-m.pos(vnn))/sqr(12*avg_hlen));
                sqlap+=lap[vnn];
                w += 1 + 1/valency(m, vnn);
                wsum += wgt;
            });
            w /= valency(m, vn);
            sqlap /= wsum;
            sqlap -= lap[vn];
            new_pos[vn] = m.pos(vn) - (0.4/w)*sqlap;
        }
        swap(new_pos, m.positions_attribute_vector());
    }
    
    return true;
}


void smooth_and_refit(HMesh::Manifold& m_in,  HMesh::Manifold& m_ref, int iter, double alpha, bool preserve_poles)
{
    auto ptav = register_faces(m_in, m_ref, preserve_poles);
    smooth_geodesic(m_in, m_ref, ptav, iter, alpha);
}


void polar_subdivide(HMesh::Manifold& mani, int MAX_ITER)
{
    
    // Time stamp all vertices
    int T=0;
    VertexAttributeVector<int> time_stamp;
    VertexAttributeVector<int> pole_tag(mani.allocated_vertices(), 0);
    for(VertexIDIterator vid = mani.vertices_begin(); vid != mani.vertices_end(); ++vid)
    {
        time_stamp[*vid] = T;
        if(is_pole(mani, *vid))
            pole_tag[*vid] = 1;
    }
    
    // For MAX_ITER subdivision steps
    for(int iter = 0; iter < MAX_ITER; ++iter)
    {
        // two iterations since we split first U edges (looping around tubes)
        // and then V edges (lying along tubes)
        for(int _uvt = 1 ; _uvt >=0; --_uvt)
        {
            HalfEdgeAttributeVector<EdgeInfo> edge_info = label_PAM_edges(mani);
            
            EdgeType et = ((_uvt == 1) ? RIB : SPINE);
            T=T+1;
            
            // Make a list of edges that are of proper type (U or V)
            vector<HalfEdgeID> split_list;
            for(HalfEdgeIDIterator hid = mani.halfedges_begin(); hid != mani.halfedges_end(); ++hid)
            {
                Walker w = mani.walker(*hid);
                if(((edge_info[*hid].is_rib() && et == RIB) ||
                    (edge_info[*hid].is_spine() && et == SPINE))
                   && w.opp().halfedge()<w.halfedge())
                    split_list.push_back(*hid);
            }
            
            // Split edges (remember to divide the UV segments.
            for(int i=0;i<split_list.size(); ++i)
                time_stamp[mani.split_edge(split_list[i])] =T;
            
            // Divide faces.
            vector<FaceID> face_list;
            for(FaceIDIterator fid=mani.faces_begin(); fid != mani.faces_end(); ++fid)
                face_list.push_back(*fid);
            for(int i=0;i<face_list.size(); ++i)
            {
                FaceID f = face_list[i];
                vector<VertexID> end_points;
                int polar_end_point = -1;
                int k=0;
                
                // Loop around a face to find the vertices that we connect when creating new faces.
                // For each face, two vertices must be selected. If a U edge has been split in a (formerly)
                // triangular face then we connect the vertex inserted on the U edge with the polar vertex.
                for(Walker w = mani.walker(f);!w.full_circle(); w = w.next())
                {
                    if(time_stamp[w.vertex()] == T)
                    {
                        end_points.push_back(w.vertex());
                        ++k;
                    }
                    else if(et == RIB && pole_tag[w.vertex()])
                    {
                        end_points.push_back(w.vertex());
                        polar_end_point = k;
                        ++k;
                    }
                }
                assert(end_points.size() == 2);
                if(k==2)
                    mani.split_face_by_edge(f, end_points[0], end_points[1]);
            }
        }
        
        // ----- CC SMOOTH
        VertexAttributeVector<Vec3d> new_vertices(mani.allocated_faces(), Vec3d(0));
        for(FaceIDIterator fi = mani.faces_begin(); fi != mani.faces_end(); ++fi)
        {
            FaceID f = *fi;
            Walker w = mani.walker(f);
            for(; !w.full_circle(); w = w.next())
            {
                VertexID v = w.vertex();
                float val = valency(mani, v);
                float A = (1.0f-3.0f/val)	* (1.0f/val);
                float B = sqr(1.0f/val);
                Walker w2 = mani.walker(f);
                for(; !w2.full_circle(); w2 = w2.next())
                {
                    VertexID v2 = w2.vertex();
                    if(v==v2)
                        new_vertices[v] += A * mani.pos(v2);
                    else
                        new_vertices[v] += B * mani.pos(v2);
                }
                
            }
        }
        
        // Polar smooth
        VertexAttributeVector<int> touched(mani.no_vertices(), 0);
        for(VertexIDIterator vi = mani.vertices_begin(); vi != mani.vertices_end(); ++vi)
        {
            if(pole_tag[*vi])
            {
                Vec3d p00 = mani.pos(*vi) * 0.5;
                Walker w = mani.walker(*vi);
                int nm = valency(mani, *vi);
                vector<Vec3d> new_p1(nm,Vec3d(0));
                for(int j=0;!w.full_circle(); w = w.circulate_vertex_ccw(),++j)
                {
                    Vec3d npos = mani.pos(w.vertex());
                    p00 += (0.5/nm)*npos;
                    
                    Walker w2 = w;
                    for(int h=0;h<nm;++h, w2=w2.circulate_vertex_ccw())
                        
                    {
                        float g = h/float(nm);
                        new_p1[j] +=
                        (1.0/nm)*(1.0 + 2.0*cos(2.0*M_PI*g)+ cos(4.0*M_PI*g)) * mani.pos(w2.vertex());
                    }
                    
                }
                
                mani.pos(*vi) = p00;
                touched[*vi] = 1;
                w = mani.walker(*vi);
                for(int j=0;!w.full_circle(); w = w.circulate_vertex_ccw(),++j)
                {
                    mani.pos(w.vertex()) = new_p1[j];
                    touched[w.vertex()] = 1;
                }
                
            }
        }
        for(VertexIDIterator vi = mani.vertices_begin(); vi != mani.vertices_end(); ++vi)
            if(!touched[*vi])
                mani.pos(*vi) = new_vertices[*vi];
        
    }
    
    
}



void polar_add_branch(HMesh::Manifold& m, HMesh::VertexAttributeVector<int>& vs)
{
    HalfEdgeID h = m.slit_edges(vs);
    FaceID f = m.close_hole(h);
    vector<Vec3d> npos;
    for(Walker w = m.walker(f); !w.full_circle(); w = w.next())
    {
        Vec3d p(0);
        p /= circulate_vertex_ccw(m, w.vertex(), [&](VertexID v) {
            p += m.pos(v);
        });
        npos.push_back(p);
    }
    int i=0;
    circulate_face_ccw(m, f, [&](VertexID v) {m.pos(v) = npos[i++];});
    m.split_face_by_vertex(f);
}

void refine_poles(HMesh::Manifold& m, HMesh::VertexAttributeVector<int>& vs)
{
    for(auto vid:m.vertices())
        if(vs[vid]==1 && is_pole(m,vid))
        {
            vector<HalfEdgeID> hvec;
            circulate_vertex_ccw(m, vid, [&](HalfEdgeID h){
                hvec.push_back(h);
            });
            for(auto h: hvec) m.split_edge(h);
            circulate_vertex_ccw(m, vid, [&](Walker w){
                m.split_face_by_edge(w.face(), w.prev().opp().vertex(), w.vertex());
            });
            
        }
}


void make_height_fun(const HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun,
                     double& vmin, double& vmax)
{
    VertexIDIterator vid = m.vertices_begin();
    VertexAttributeVector<int> nailed(m.allocated_vertices(), 0);
    vmin = vmax = m.pos(*vid)[1];
    for(; vid != m.vertices_end(); ++vid)
    {
        double v = dot(m.pos(*vid),Vec3d(0.0,1,0.00));
        fun[*vid] = v;
        vmin = min(v, vmin);
        vmax = max(v, vmax);
    }
    
}

void make_adf_fun(HMesh::Manifold& m, double t, HMesh::VertexAttributeVector<double>& F,
                  double& vmin, double& vmax)
{
    static Harmonics harm(m);
    harm.compute_adf(F, t, 0);
    VertexAttributeVector<int> nailed(m.allocated_vertices(), 0);
    vmin = vmax = F[*m.vertices_begin()];
    for(auto vid : m.vertices())
    {
        vmin = min(F[vid], vmin);
        vmax = max(F[vid], vmax);
    }
    cout << vmin << " " << vmax << endl;
}

