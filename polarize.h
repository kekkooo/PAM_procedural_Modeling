//
//  polarize.h
//  GEL
//
//  Created by J. Andreas Bærentzen on 18/03/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <GEL/CGLA/Vec2d.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/AttributeVector.h>

#ifndef __POLARIZE_H__
#define __POLARIZE_H__

typedef std::pair<double, std::vector<HMesh::HalfEdgeID>> Region;
enum EdgeType
{
    UNKNOWN, SPINE, RIB, RIB_JUNCTION
};

struct EdgeInfo {
    EdgeType edge_type;
    int id;
    
    bool is_rib()       const { return edge_type == RIB || edge_type == RIB_JUNCTION;}
    bool is_spine()     const { return edge_type == SPINE;}
    bool is_junction()  const { return edge_type == RIB_JUNCTION; }
    
    EdgeInfo(): edge_type(UNKNOWN), id(0) {}
    EdgeInfo(EdgeType _edge_type, int _id): edge_type(_edge_type), id(_id) {}
};

HMesh::HalfEdgeAttributeVector<EdgeInfo> label_PAM_edges(HMesh::Manifold& m);
	


const CGLA::Vec3f& get_color(int i);

bool is_pole(HMesh::Manifold& m, HMesh::VertexID v);

void make_height_fun(const HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun,
                     double& vmin, double &vmax);

void make_adf_fun(HMesh::Manifold& m, double t, HMesh::VertexAttributeVector<double>& fun,
                  double& vmin, double& vmax);

template<typename T>
void smooth_fun(const HMesh::Manifold& m,
                const HMesh::VertexAttributeVector<int>& nailed,
                HMesh::VertexAttributeVector<T>& fun, int iter);

void polarize_mesh(HMesh::Manifold& m, HMesh::VertexAttributeVector<double>& fun, double vmin, double vmax, const int divisions);

void show_skin(HMesh::Manifold& m, int j);
void skeleton_retract(HMesh::Manifold& m, double frac, HMesh::VertexID v=HMesh::InvalidVertexID);

void polar_subdivide(HMesh::Manifold& mani, int MAX_ITER);
void polar_segment(HMesh::Manifold& m, bool show_segments=false);
int number_rib_edges(HMesh::Manifold& m,  HMesh::HalfEdgeAttributeVector<EdgeInfo>& edge_info, HMesh::HalfEdgeID h);


Region trace_region(HMesh::Manifold& m, HMesh::HalfEdgeID hid);
void collapse_region(HMesh::Manifold& m, Region& region);
void refine_region(HMesh::Manifold& m, Region& region);

void simplify_polar_mesh(HMesh::Manifold& m, double frac, int max_iter=1);
void smooth_and_refit(HMesh::Manifold& m1,  HMesh::Manifold& m2, int iter, double alpha, bool preserve_poles);
void polar_add_branch(HMesh::Manifold& m, HMesh::VertexAttributeVector<int>& vs);
void refine_poles(HMesh::Manifold& m, HMesh::VertexAttributeVector<int>& vs);

void dist_color_vertices(HMesh::Manifold& m, HMesh::VertexAttributeVector<int>& vs);
bool pole_quadrangulation(HMesh::Manifold& m);

void skeletonize(HMesh::Manifold& m, const HMesh::VertexAttributeVector<double>& fun, double w=1.0);
void skeleton_line_field(HMesh::Manifold& m, const HMesh::VertexAttributeVector<double>& fun, HMesh::VertexAttributeVector<CGLA::Vec3d>& lin);
#endif
