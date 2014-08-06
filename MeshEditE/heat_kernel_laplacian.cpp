//
//  heat_kernel_laplacian.cpp
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 05/02/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "LogMap.h"
#include "heat_kernel_laplacian.h"
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include <GEL/LinAlg/LapackFunc.h>

using namespace CGLA;
using namespace HMesh;
using namespace Geometry;
using namespace GLGraphics;
using namespace LinAlg;
using namespace std;

VertexAttributeVector<double> heat_kernel_laplacian(Manifold& m, double R) {
    double avg_len = 0;
    for(auto h : m.halfedges())
        avg_len += length(m,h);
    avg_len /= m.no_halfedges();
    R *= avg_len;
    
    int cnt = 0;
    double avg = 0;
    VertexAttributeVector<double> mc(m.allocated_vertices(), 0);
    for(VertexID v: m.vertices()) {
        FaceID f = m.walker(v).face();
        Vec3f b = barycentric_coords(m, m.pos(v), f);
        LogMap log_map(m, ManiPoint(f,b), R);
        auto pts = log_map.enumerate();
        int N = pts.size();
        
        vector<double> weight_vec(N);   
        CMatrix U(N,5);
        for(int i=0;i<N;++i)
        {
            Vec2f uv = pts[i].second;
            U[i][0] = uv[0];
            U[i][1] = uv[1];
            U[i][2] = 0.5*sqr(uv[0]);
            U[i][3] = uv[0]*uv[1];
            U[i][4] = 0.5*sqr(uv[1]);
        }
        
        CMatrix UT = U.Transposed();
        CMatrix D = Inverted(UT * U)*UT;
        for(int i=0;i<N;++i)
            weight_vec[i] = D[2][i] + D[4][i];

        
        Vec3d p03d = m.pos(v);
        Vec3d L(0);
        for(int i=0;i<N;++i) {
            Vec3d p3d = pts[i].first;
            L += weight_vec[i] * (p3d-p03d);
        }
        mc[v] = - 0.5 * L.length() * sign(dot(L,normal(m,v)));
        avg += mc[v];
        ++cnt;
    }
    cout << "Average mean curvature : " << avg/cnt << endl;
    return mc;
}

//VertexAttributeVector<double> heat_kernel_laplacian(Manifold& m, int N, double R) {
//    double avg_len = 0;
//    for(auto h : m.halfedges())
//        avg_len += length(m,h);
//    avg_len /= m.no_halfedges();
//    R *= avg_len;
//    
//    vector<Vec2f> pts(N);
//    vector<double> weight_vec(N);
//    
//    for(int i=0;i<N;++i)
//    {
//        Vec2f p;
//        do {
//            p = Vec2f(gel_rand(), gel_rand());
//            p /= GEL_RAND_MAX;
//            p *= 2.0;
//            p -= Vec2f(1);
//        } while(sqr_length(p)>1.0);
//        pts[i] = R*p;
//        
//    }
//    CMatrix U(N,5);
//    for(int i=0;i<N;++i)
//    {
//        U[i][0] = pts[i][0];
//        U[i][1] = pts[i][1];
//        U[i][2] = 0.5*sqr(pts[i][0]);
//        U[i][3] = pts[i][0]*pts[i][1];
//        U[i][4] = 0.5*sqr(pts[i][1]);
//    }
//    
//    CMatrix UT = U.Transposed();
//    CMatrix D = Inverted(UT * U)*UT;
//    for(int i=0;i<N;++i)
//        weight_vec[i] = D[2][i] + D[4][i];
//    
//    int cnt = 0;
//    double avg = 0;
//    VertexAttributeVector<double> mc(m.allocated_vertices(), 0);
//    for(VertexID v: m.vertices()) {
//        FaceID f = m.walker(v).face();
//        Vec3f b = barycentric_coords(m, m.pos(v), f);
//        LogMap log_map(m, ManiPoint(f,b), 3 * R);
//        Vec3d p03d = m.pos(v);
//        Vec3d L(0);
//        for(int i=0;i<N;++i) {
//            Vec3d p3d = manipoint_to_3D(m, log_map.find_uv(pts[i]));
//            L += weight_vec[i] * (p3d-p03d);
//        }
//        mc[v] = - 0.5 * L.length() * sign(dot(L,normal(m,v)));
//        avg += mc[v];
//        ++cnt;
//    }
//    cout << "Average mean curvature : " << avg/cnt << endl;
//    return mc;
//}
//
//
//
