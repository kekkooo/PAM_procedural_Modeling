//
//  LogMap.cpp
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 06/01/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include <thread>
#include <queue>
#include <GEL/HMesh/curvature.h>
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include "LogMap.h"
#include "HMeshParallelKit.h"
#include <GEL/Geometry/KDTree.h>

using namespace CGLA;
using namespace Geometry;
using namespace GLGraphics;
using namespace std;
using namespace HMesh;

ManiPoint LogMap::find_uv(const Vec2f& uv) {
    ManiPoint pt;
    FaceAttributeVector<int> visited(m.allocated_faces(), 0);
    queue<FaceID> Q;
    Q.push(f0);
    while(!Q.empty())
    {
        pt.f = Q.front();
        visited[pt.f] = 1;
        Q.pop();
        pt.b = uv_barycentric(pt.f, uv);
        if(pt.b[0]>=0 && pt.b[1]>=0 && pt.b[2]>=0)
            return pt;
        int i=1;
        circulate_face_ccw(m, pt.f, [&](FaceID fn) {
            if(!visited[fn] && pt.b[i]<=0)
                Q.push(fn);
            i = (i+1)%3;
        });
    }
    cout << "FAILED for " << uv <<  endl;
    pt.b = Vec3f(1.0/3);
    pt.f = f0;
    return pt;
}

vector<pair<Vec3d, Vec2f>> LogMap::enumerate() {
    vector<pair<Vec3d, Vec2f>> pts;
    VertexAttributeVector<int> visited(m.allocated_faces(), 0);
    queue<VertexID> Q;
    Q.push(m.walker(f0).vertex());
    while(!Q.empty())
    {
        VertexID v = Q.front();
        Q.pop();
        visited[v] = 1;
        pts.push_back(make_pair(m.pos(v), uv_coords(v)));
        circulate_vertex_ccw(m, v, [&](VertexID vn) {
            if(polar_map[vn][0]<DBL_MAX && !visited[vn])
                Q.push(vn);
        });
    }
    return pts;
}



void LogMap::init_vertex(VertexID s)
{
    double angle_sum = 0;
    Vec3d s_pos = m.pos(s);
    circulate_vertex_ccw(m, s, [&](Walker w){
        VertexID v = w.vertex();
        Vec3d e0 = m.pos(v) - s_pos;
        polar_map[v][0] = length(e0);
        Vec3d e1 = m.pos(w.prev().opp().vertex()) - s_pos;
        e0.cond_normalize();
        e1.cond_normalize();
        polar_map[v][1] = angle_sum;
        double a = acos(min(1.0, max(-1.0, dot(e0,e1))));
        angle_sum += a;
    });
    
    circulate_vertex_ccw(m, s, [&](VertexID v){
        polar_map[v][1] *= 2.0 * M_PI / angle_sum;
    });
}

void LogMap::init_face(FaceID f, const Vec3f& bary)
{
    Vec3d s_pos(0);
    int i =0;
    circulate_face_ccw(m, f, [&](VertexID v) { s_pos += bary[i++]*m.pos(v);});
    
    double angle_sum = 0;
    circulate_face_ccw(m, f, [&](Walker w) {
        VertexID v = w.vertex();
        Vec3d dir = m.pos(v)-s_pos;
        polar_map[v][0] = length(dir);
        polar_map[v][1] = angle_sum;
        Vec3d dirn = m.pos(w.next().vertex()) - s_pos;
        angle_sum += acos(min(1.0,max(-1.0, dot(dir, dirn)/(dir.length()*dirn.length()))));
    });
}

double LogMap::compute_distance(VertexID s, HalfEdgeID h, double& alpha)
{
    Vec3d pt = m.pos(s);
    Walker w = m.walker(h);
    VertexID v0 = w.opp().vertex();
    VertexID v1 = w.vertex();
    const Vec3d& Nk = m.pos(v1);
    const Vec3d& Nj = m.pos(v0);
    
    const double Uk = polar_map[v1][0];
    const double Uj = polar_map[v0][0];
    
    const double djptsq = sqr_length(Nj-pt);
    const double dkptsq = sqr_length(Nk-pt);
    
    const Vec3d ekj = Nk-Nj;
    const double djksq = sqr_length(ekj);
    const double djk = sqrt(djksq);
    
    //Stable evaluation of Herons formula
    const double a = max(djk, max(Uj, Uk));
    const double c = min(djk, min(Uj, Uk));
    double b;
    if(a == djk)
        b = max(Uj, Uk);
    else if(a == Uj)
        b = max(djk, Uk);
    else
        b = max(djk, Uj);
    
    const double H_under_root = ( (a + (b+c)) *
                                 (c - (a-b)) *
                                 (c + (a-b)) *
                                 (a + (b-c)) );
    
    if(H_under_root < 0 || djk < 10e-12) {
        // Triangle inequality fails, return Dijkstra instead
        const double dijkstra_j = Uj + length(Nj-pt);
        const double dijkstra_k = Uk + length(Nk-pt);
        if(dijkstra_j < dijkstra_k) {
            alpha = 0;
            return dijkstra_j;
        } else {
            alpha = 1;
            return dijkstra_k;
        }
    }
    
    
    const double H = sqrt( H_under_root );
    
    const Vec3d ej = Nj-pt;
    const Vec3d ek = Nk-pt;
    
    const double A2 = length(cross(ej,ek));
    
    const double ej_ekj = dot(ej, ekj);
    const double ek_ekj = dot(ek, ekj);
        
    const double Ujsq = Uj*Uj;
    const double Uksq = Uk*Uk;
    
    const double f1_j = A2 * (djksq + Uksq - Ujsq);
    const double f1_k = A2 * (djksq + Ujsq - Uksq);
    
    const double xj = (f1_j + ek_ekj*H);
    const double xk = (f1_k - ej_ekj*H);
    
    if(xj < 0 || xk < 0) {
        // Update from outside triangle, return Dijkstra instead
        const double dijkstra_j = Uj + length(Nj-pt);
        const double dijkstra_k = Uk + length(Nk-pt);
        if(dijkstra_j < dijkstra_k) {
            alpha = 0;
            return dijkstra_j;
        } else {
            alpha = 1;
            return dijkstra_k;
        }
    }
    
    const double f4 = 2*A2*djksq;
    
    double Ui = sqrt(xj*xj*djptsq + 2*xj*xk*dot(ej,ek) + xk*xk*dkptsq)/f4;
    
    const double cos_jk = (  Ujsq+Uksq - djksq)/(2*Uj*Uk);
    const double cos_ji  = ( Ujsq+Ui*Ui - djptsq)/(2*Uj*Ui);
    alpha = acos(cos_ji) /  acos (cos_jk);
    
    return Ui;
}

double LogMap::compute_angle(VertexID s, HalfEdgeID h, double alpha)
{
    Walker w = m.walker(h);
    double nkphi = polar_map[w.vertex()][1];
    double njphi = polar_map[w.opp().vertex()][1];
    
    const double diff = fabs(njphi-nkphi);
    
    if(diff < 1e-6)
        return njphi;
    
    if(diff > M_PI) {
        //Make the interpolation modulo 2pi
        if(njphi < nkphi) njphi += 2*M_PI;
        else nkphi += 2*M_PI;
    }
    
    double angle = (1-alpha) * njphi  + alpha * nkphi;
    
    if(angle > 2*M_PI)
        angle -= 2*M_PI;
    
    return angle;
}

LogMap::LogMap(const Manifold& _m, ManiPoint pt, double max_dist): m(_m), f0(pt.f)
{
    polar_map = VertexAttributeVector<Vec2d>(m.allocated_vertices(), Vec2d(DBL_MAX,-1));
    priority_queue<pair<double,VertexID>> Q;
    
    
    if(pt.b.max_coord()>0.999) {
        VertexID s;
        int i=0;
        circulate_face_ccw(m, f0, [&](VertexID v){ if(pt.b[i++]>0.999) s = v; });
        init_vertex(s);
        polar_map[s][0] = 0.0;
        polar_map[s][1] = 0.0;
        circulate_vertex_ccw(m, s, [&](VertexID v) {Q.push(make_pair(-polar_map[v][0], v));});
    }
    else {
        if(pt.b.min_coord()<0.001) {
            pt.b += Vec3f(0.1);
            pt.b /= pt.b[0]+pt.b[1]+pt.b[2];
        }
        init_face(f0, pt.b);
        circulate_face_ccw(m, f0, [&](VertexID v) {Q.push(make_pair(-polar_map[v][0], v));});
    }
    
    while(!Q.empty())
    {
        VertexID j = Q.top().second;
        Q.pop();
        
        circulate_vertex_ccw(m, j, [&](VertexID i){
            // Compute distance
            double alpha;
            HalfEdgeID he_id;
            double new_dist_i=DBL_MAX;
            circulate_vertex_ccw(m, i, [&](Walker w){
                w = w.next();
                if(polar_map[w.vertex()][0] < DBL_MAX &&
                   polar_map[w.opp().vertex()][0] < DBL_MAX) {
                    double a;
                    HalfEdgeID h = w.halfedge();
                    double d = compute_distance(i, h , a);
                    if(d<new_dist_i) {
                        new_dist_i = d;
                        alpha = a;
                        he_id = h;
                    }
                }
            });
            if(polar_map[i][0]/new_dist_i > (1.000000001))
            {
                polar_map[i][0] = new_dist_i;
                polar_map[i][1] = compute_angle(i, he_id, alpha);
                if(new_dist_i < max_dist)
                    Q.push(make_pair(-new_dist_i, i));
            }
        });
    }
}

void dist_color_vertices(Manifold& m, HMesh::VertexAttributeVector<int>& vs)
{
    double avglen = 0;
    for(HalfEdgeID h: m.halfedges()) avglen += length(m,h);
    avglen /= m.no_halfedges();
    
    for(VertexID v: m.vertices())
        if(vs[v]) {
            Walker w = m.walker(v);
            ManiPoint pt(w.face(), Vec3f(1,0,0));
            LogMap log_map(m, pt, 10.0*avglen);
            for(VertexID v: m.vertices())
                CheckerBoardRenderer::param[v] = log_map.uv_coords(v) / avglen;
            break;
        }
    
}

inline bool finite(const CGLA::Vec2f& v) {
    return v[0] <= FLT_MAX && v[0] >= -FLT_MAX &&
    v[1] <= FLT_MAX && v[1] >= -FLT_MAX;
}


Vec3f barycentric_coords(Manifold& m, const Vec3d& _p, FaceID f) {
    Vec3d c = centre(m,f);
    Vec3d n = normal(m, f);
    Vec3d p = _p-n*dot(n,_p-c);
    
    Vec3f bary(0);
    int i = 0;
    circulate_face_ccw(m, f, [&](Walker w) {
        Vec3d uv1 = m.pos(w.prev().vertex());
        Vec3d uv0 = m.pos(w.prev().opp().vertex());
        bary[i++] = length(cross(uv1-uv0,p-uv0));
    });
    bary /= dot(Vec3f(1), bary);
    return bary;
}

Vec3d manipoint_to_3D(Manifold& m, const ManiPoint& pt)
{
    Walker w = m.walker(pt.f);
    array<Vec3d,3> p = {m.pos(w.vertex()),
        m.pos(w.next().vertex()),
        m.pos(w.next().next().vertex())};
    Vec3d q(0);
    for(int i : {0,1,2})
        q += pt.b[i] * p[i];
    return q;
//    array<Vec3d,3> n = {normal(m,w.vertex()),
//        normal(m, w.next().vertex()),
//        normal(m, w.next().next().vertex())};
//    Vec3d q2(0);
//    for(int i : {0,1,2})
//        q2 += pt.b[i] * (q-n[i]*dot(n[i],q-p[i]));
//    return q2;
}

VertexAttributeVector<ManiPoint> register_faces(HMesh::Manifold& m_in,  HMesh::Manifold& m_ref,
                                                bool preserve_poles)
{
    auto is_pole = [](Manifold& m, VertexID v) -> bool
    {
        Walker w = m.walker(v);
        for(; !w.full_circle(); w = w.circulate_vertex_ccw())
            if(w.face() == InvalidFaceID || no_edges(m, w.face()) != 3)
                return false;
        return true;
    };

    KDTree<Vec3d, FaceID> face_cloud;
    
    double avgsz =0;
    for(auto fid : m_ref.faces())
        avgsz += area(m_ref, fid);
    avgsz /= m_ref.no_faces();
    
    for(auto fid : m_ref.faces()) {
        face_cloud.insert(centre(m_ref,fid), fid);
        double A = area(m_ref, fid);
        for(int i=0;i<int(A/avgsz);++i)
        {
            Vec3d b = Vec3d(gel_rand(), gel_rand(),gel_rand());
            b /= dot(Vec3d(1), b);
            Vec3d p(0);
            int j=0;
            circulate_face_ccw(m_ref, fid, [&](VertexID v){
                p += b[j++] * m_ref.pos(v);
            });
            face_cloud.insert(p, fid);
        }
    }
    face_cloud.build();
    
    Vec3d c;
    float r;
    bsphere(m_ref, c, r);
    r *= 0.3;
    
    VertexAttributeVector<ManiPoint> ptav(m_in.allocated_vertices());
    for(auto vid : m_in.vertices()) {
        double dist = r;
        Vec3d p = m_in.pos(vid);
        Vec3d key;
        FaceID f = InvalidFaceID;
        
        face_cloud.closest_point(p, dist, key, f);
        Vec3f bary = barycentric_coords(m_ref, p, f);
        if(preserve_poles && is_pole(m_in, vid))
            ptav[vid] = ManiPoint(f, bary, true);
        else
            ptav[vid] = ManiPoint(f,bary);
    }
    return ptav;
}

void manipoints_to_3D(Manifold& m_in, Manifold& m_ref, VertexAttributeVector<ManiPoint>& ptav) {
    for(VertexID v: m_in.vertices())
        if(!ptav[v].fixed)
            m_in.pos(v) = manipoint_to_3D(m_ref, ptav[v]);
}


void smooth_geodesic(Manifold& m_in, Manifold& m_ref, VertexAttributeVector<ManiPoint>& ptav, int max_iter, float weight)
{
    auto vertex_ids = batch_vertices(m_in);
    
    VertexAttributeVector<ManiPoint> ptav_new = ptav;
    
    auto f = [&](const vector<VertexID>& vids) {
        for(VertexID v: vids)
            if(!ptav[v].fixed)
            {
                double len=0;
                circulate_vertex_ccw(m_in, v, [&](HalfEdgeID h){len=max(len,length(m_in,h));});
                LogMap log_map(m_ref, ptav[v], 3.5 * len);
                Vec2f uv = log_map.barycentric_uv(ptav[v]);
                Vec2f new_uv(0);
                //int cnt=0;
                double wsum =0;
                circulate_vertex_ccw(m_in, v, [&](VertexID vn){
                    Vec2f uvn = log_map.barycentric_uv(ptav[vn]);
                    if(finite(uvn)) {
                        double w = 1;//sqr_length(uv-uvn)/sqr_length(m_in.pos(v)-m_in.pos(vn));
                        wsum += w;
                        new_uv += w*uvn;
//                        ++cnt;
                    }
                });
                new_uv = (1-weight) * uv + weight * new_uv/wsum;
                ptav_new[v] = log_map.find_uv(new_uv);
            }
    };
    
    for(auto _ : range(0, max_iter))  {
        for_each_vertex_parallel(CORES, vertex_ids, f);
        swap(ptav, ptav_new);
        cout << "." << flush;
    }
    manipoints_to_3D(m_in, m_ref, ptav);
}

void smooth_geodesic(Manifold& m_in, Manifold& m_ref, int iter, float weight)
{
    VertexAttributeVector<ManiPoint> ptav = register_faces(m_in, m_ref);
    
//    for(HalfEdgeID h: m_in.halfedges())
//        if(m_in.in_use(h))
//    {
//        double len=0;
//        Walker w =m_in.walker(h);
//        
//        VertexID v = w.vertex();
//        circulate_vertex_ccw(m_in, v, [&](HalfEdgeID h){len=max(len,length(m_in,h));});
//        LogMap log_map(m_ref, ptav[v], 3.5 * len);
//        
//        Vec2f uv = log_map.barycentric_uv(ptav[v]);
//        Vec2f uvo =log_map.barycentric_uv(ptav[w.opp().vertex()]);
//        
//        Vec2f uvn = log_map.barycentric_uv(ptav[w.next().vertex()]);
//        Vec2f uvon = log_map.barycentric_uv(ptav[w.opp().next().vertex()]);
//        
//        Vec3d p = m_in.pos(v);
//        Vec3d po = m_in.pos(w.opp().vertex());
//        Vec3d pn = m_in.pos(w.next().vertex());
//        Vec3d pon = m_in.pos(w.opp().next().vertex());
//        
//        double R = sqr_length(uv-uvo)/sqr_length(p-po);
//        double Rf = sqr_length(uvn-uvon)/sqr_length(pn-pon);
////        if(Rf<0.8*R && h<w.opp().halfedge() && precond_flip_edge(m_in, h))
////            m_in.flip_edge(h);
////        else
//        {
//            Vec2f uvc = (uv + uvo + uvn)/3.0;
//            ManiPoint mp = log_map.find_uv(uvc);
//            Vec3d pc = manipoint_to_3D(m_ref, mp);
//            Vec3d n = normalize(cross(pn-p,po-p));
//            double d = abs(dot(normalize(pc-p), n));
//            if(d>0.4)
//            {
//                cout << "splitting" << endl;
//                VertexID vnew = m_in.split_face_by_vertex(w.face());
//                m_in.pos(vnew) = pc;
//                ptav[vnew] = mp;
//            }
//            
//        }
//    }
    
    
    
    //    for(VertexID v: m_in.vertices()) {
    //        double a = 0;
    //        circulate_vertex_ccw(m_in, v, [&](FaceID f){a+=area(m_in, f);});
    //        double g = gaussian_curvature_angle_defect(m_in, v)*a;
    //
    //        if(g<-M_PI/16)
    //            ptav[v].fixed = true;
    //    }
    smooth_geodesic(m_in,m_ref,ptav, iter, weight);
}
