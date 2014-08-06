//
//  LogMap.h
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 06/01/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__LogMap__
#define __MeshEditE__LogMap__

#include <iostream>
#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/Vec2f.h>
#include <GEL/CGLA/Vec3f.h>
#include <GEL/HMesh/Manifold.h>

struct ManiPoint
{
    HMesh::FaceID f;
    CGLA::Vec3f b;
    bool fixed;
    
    ManiPoint(): f(HMesh::InvalidFaceID), b(0), fixed(false) {}
    ManiPoint(HMesh::FaceID _f, const CGLA::Vec3f& _b, bool _fixed = false): f(_f), b(_b), fixed(_fixed) {}
    bool good() {return f != HMesh::InvalidFaceID;}
};


class LogMap
{
    const HMesh::Manifold& m;
    HMesh::VertexAttributeVector<CGLA::Vec2d> polar_map;
    HMesh::FaceID f0;
    
    void init_vertex(HMesh::VertexID s);
    void init_face(HMesh::FaceID f, const CGLA::Vec3f& bary);
    double compute_distance(HMesh::VertexID s, HMesh::HalfEdgeID h, double&  alpha);
    double compute_angle(HMesh::VertexID s, HMesh::HalfEdgeID h, double alpha);
    
    
public:
    LogMap(const HMesh::Manifold& m, ManiPoint pt, double max_dist);
    CGLA::Vec2f uv_coords(HMesh::VertexID v) {
        return polar_map[v][0] * CGLA::Vec2f(cos(polar_map[v][1]),sin(polar_map[v][1]));
    }
    
    CGLA::Vec3f uv_barycentric(HMesh::FaceID f, const CGLA::Vec2f& uv)
    {
        CGLA::Vec3f bary(0);
        int i = 0;
        circulate_face_ccw(m, f, [&](HMesh::Walker w) {
            CGLA::Vec2f uv1 = uv_coords(w.prev().vertex());
            CGLA::Vec2f uv0 = uv_coords(w.prev().opp().vertex());
            bary[i++] = cross(uv1-uv0,uv-uv0);
        });
        return bary/dot(bary,CGLA::Vec3f(1));
    }
    
    CGLA::Vec2f barycentric_uv(const ManiPoint& pt)
    {
        HMesh::Walker w = m.walker(pt.f);
        return pt.b[0] * uv_coords(w.vertex()) +
        pt.b[1] * uv_coords(w.next().vertex()) +
        pt.b[2] * uv_coords(w.next().next().vertex());
    }
    
    ManiPoint find_uv(const CGLA::Vec2f& uv);
    
    std::vector<std::pair<CGLA::Vec3d, CGLA::Vec2f>> enumerate();
};

CGLA::Vec3f barycentric_coords(HMesh::Manifold& m, const CGLA::Vec3d& _p, HMesh::FaceID f);
CGLA::Vec3d manipoint_to_3D(HMesh::Manifold& m, const ManiPoint& pt);
void manipoints_to_3D(HMesh::Manifold& m_in, HMesh::Manifold& m_ref,
                      HMesh::VertexAttributeVector<ManiPoint>& ptav);

void dist_color_vertices(HMesh::Manifold& m, HMesh::VertexAttributeVector<int>& vs);

HMesh::VertexAttributeVector<ManiPoint> register_faces(HMesh::Manifold& m_in,  HMesh::Manifold& m_ref,
                                                bool preserve_poles = false);

void smooth_geodesic(HMesh::Manifold& m_in, HMesh::Manifold& m_ref, HMesh::VertexAttributeVector<ManiPoint>& ptav, int max_iter=20, float weight=0.5);

void smooth_geodesic(HMesh::Manifold& m_in, HMesh::Manifold& m_ref, int max_iter=20, float weight=0.5);



#endif /* defined(__MeshEditE__LogMap__) */
