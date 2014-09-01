//
//  additional_console_funcs.cpp
//  MeshEditE
//
//  Created by J. Andreas Bærentzen on 14/10/13.
//  Copyright (c) 2013 J. Andreas Bærentzen. All rights reserved.
//

#include <GEL/HMesh/Manifold.h>
#include <GEL/HMesh/smooth.h>
#include <GEL/HMesh/AttributeVector.h>
#include <GEL/HMesh/mesh_optimization.h>
#include <GEL/HMesh/triangulate.h>
#include <strstream>
#include <istream>
#include <fstream>
#include <string>
#include <regex>
#include <set>
#include <queue>
#include "additional_console_funcs.h"
#include <GEL/GLGraphics/MeshEditor.h>

#include "LogMap.h"
#include "polarize.h"
#include "heat_kernel_laplacian.h"

using namespace std;
using namespace CGLA;
using namespace GLGraphics;
using namespace HMesh;

void root3_subdivide(Manifold& m_in, Manifold& m, VertexAttributeVector<double>& fun)
{
    if(&m != &m_in)
        m = m_in;
    
    VertexAttributeVector<int> vtouched(m.allocated_vertices(), 0);
    VertexAttributeVector<Vec3d> new_pos(m.allocated_vertices(), Vec3d(0));
    
    for (VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid) {
        int v = valency(m, *vid);
        double beta = (4.0 - 2.0 * cos(2.0 * M_PI / v))/(9.0*v);
        new_pos[*vid] = (1.0 - v * beta) * m.pos(*vid);
        for(Walker w = m.walker(*vid); !w.full_circle(); w = w.circulate_vertex_ccw())
        {
            new_pos[*vid] += beta * m.pos(w.vertex());
        }
    }
    
    vector<FaceID> faces;
    for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f)
        faces.push_back(*f);
    for(int i=0;i<faces.size(); ++i) {
        double f_val = 0;
        f_val /= circulate_face_ccw(m, faces[i], [&](VertexID v){f_val += fun[v];});
        VertexID vid = m.split_face_by_vertex(faces[i]);
        vtouched[vid] = 1;
        fun[vid] = f_val;
    }
    for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end(); ++h)
    {
        Walker w = m.walker(*h);
        
        if(vtouched[w.vertex()] == 0 && vtouched[w.opp().vertex()] == 0 &&
           precond_flip_edge(m, *h))
            m.flip_edge(*h);
    }
    
    for (VertexIDIterator vid = m.vertices_begin(); vid != m.vertices_end(); ++vid)
        if(vtouched[*vid] == 0)
            m.pos(*vid) = new_pos[*vid];
}



void console_simplify_polar(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    
    me->save_active_mesh();
    simplify_polar_mesh(me->active_mesh(), 0, iter);
}

void console_polar_segment(MeshEditor* me, const std::vector<std::string> & args)
{
    int segments = 1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> segments;
    }
    polar_segment(me->active_mesh(), segments);
}

void console_polar_skeleton(MeshEditor* me, const std::vector<std::string> & args)
{
    double frac = 0.9;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> frac;
    }
    VertexID v0 = InvalidVertexID;
    for(VertexID v: me->active_mesh().vertices())
        if(me->get_vertex_selection()[v]==1)
        {
            v0 = v;
            break;
        }
    
    me->save_active_mesh();
    skeleton_retract(me->active_mesh(), frac, v0);
}

void console_polar_show_skin(MeshEditor* me, const std::vector<std::string> & args)
{
    int j=0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> j;
    }
    show_skin(me->active_mesh(), j);
}


void console_polar_collapse_selected(MeshEditor* me, const std::vector<std::string> & args)
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    auto selection = me->get_vertex_selection();
    for(HalfEdgeID h: m.halfedges()) {
        Walker w = m.walker(h);
        if(h<w.opp().halfedge() &&
           selection[w.vertex()] == 1 &&
           selection[w.opp().vertex()] == 1) {
            Region r = trace_region(m, h);
            collapse_region(m, r);
        }
    }
}
void console_polar_refine_selected(MeshEditor* me, const std::vector<std::string> & args)
{
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    auto selection = me->get_vertex_selection();
    for(HalfEdgeID h: m.halfedges()) {
        Walker w = m.walker(h);
        if(h<w.opp().halfedge() &&
           selection[w.vertex()] == 1 &&
           selection[w.opp().vertex()] == 1) {
            Region r = trace_region(m, h);
            refine_region(m, r);
        }
    }
}


void console_polar_subdivide(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    me->save_active_mesh();
    polar_subdivide (me->active_mesh(), iter);
}

void console_polar_add_branch(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    me->save_active_mesh();
    polar_add_branch(me->active_mesh(), me->get_vertex_selection());
}

void console_refine_poles(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    me->save_active_mesh();
    refine_poles(me->active_mesh(), me->get_vertex_selection());
}

void bridge_one_rings(Manifold& mani, VertexID vid0, VertexID vid1)
{
	vector<VertexID> loop0;
	Walker hw0 = mani.walker(vid0);
	for(;!hw0.full_circle(); hw0 = hw0.circulate_vertex_ccw()) {
		loop0.push_back(hw0.vertex());
	}
    
	vector<VertexID> loop1;
	Walker hw1 = mani.walker(vid1);
	for(;!hw1.full_circle(); hw1 = hw1.circulate_vertex_ccw()){
		loop1.push_back(hw1.vertex());
    }
    
    vector<pair<VertexID, VertexID> > connections;
    
	size_t L0= loop0.size();
	size_t L1= loop1.size();
    
    assert(L0==L1);
    
    size_t L = L0;
    
    float min_len = FLT_MAX;
    int j_off_min_len = -1;
    for(int j_off = 0; j_off < L; ++j_off)
    {
        float len = 0;
        for(int i=0;i<L;++i)
            len += sqr_length(mani.pos(loop0[i]) - mani.pos(loop1[(L+j_off - i)%L]));
        if(len < min_len)
        {
            j_off_min_len = j_off;
            min_len = len;
        }
    }
    for(int i=0;i<L;++i)
        connections.push_back(pair<VertexID, VertexID>(loop0[i],loop1[(L+ j_off_min_len - i)%L]));
	// Merge the two one rings producing two faces.
	FaceID f0 = mani.merge_one_ring(vid0);
	FaceID f1 = mani.merge_one_ring(vid1);
	
	// Bridge the just created faces.
	vector<HalfEdgeID> newhalfedges = mani.bridge_faces(f0, f1, connections);
}


void console_merge_poles(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    auto selection = me->get_vertex_selection();
    
    int n = 0;
    VertexID v0, v1;
    for(VertexID v : m.vertices())
    {
        if(selection[v])
        {
            if(n==0)
            {
                v0 = v;
                ++n;
            }
            else{
                v1 = v;
                break;
            }
        }
    }
    if(valency(m, v0) == valency(m, v1))
        bridge_one_rings(m, v0, v1);
    else {
        me->printf("Mismatched valencies");
    }
}


void console_dist_color(MeshEditor* me, const std::vector<std::string> & args)
{
    dist_color_vertices(me->active_mesh(), me->get_vertex_selection());
}

typedef set<VertexID> VertexSet;
typedef set<FaceID> FaceSet;
typedef set<HalfEdgeID> HalfEdgeSet;

bool remove_wedge(Manifold& m, VertexSet& vset)
{
    vector<pair<int, vector<VertexID>>> path_vec;
    for(VertexID v:vset)
    {
        if(valency(m, v) == 3)
        {
            circulate_vertex_ccw(m, v, [&](Walker w){
                vector<VertexID> path;
                path.push_back(v);
                while(vset.count(v) && valency(m, w.vertex()) == 4)
                {
                    path.push_back(w.vertex());
                    w = w.next().opp().next();
                }
                if(vset.count(v) && valency(m, w.vertex())==3) {
                    path.push_back(w.vertex());
                    
                    bool path_is_good=true;
                    for(auto v: path)
                        circulate_vertex_ccw(m, v, [&](VertexID vn){
                            if(vset.count(vn)==0)
                                path_is_good=false;
                        });
                    if(path_is_good)
                        path_vec.push_back(make_pair(path.size(), path));
                }
            });
        }
    }
    if(!path_vec.empty()) {
        Manifold m_safe = m;

        sort(begin(path_vec), end(path_vec));
        
        VertexAttributeVector<int> sel_set(0);
        for(auto v: path_vec.front().second)
            sel_set[v] = 1;
        vector<HalfEdgeID> hevec;
        for(HalfEdgeID h: m.halfedges())
            if(m.in_use(h))
            {
                Walker w = m.walker(h);
                if(sel_set[w.vertex()] && !sel_set[w.opp().vertex()])
                    hevec.push_back(h);
            }
        
        vector<VertexID> garbage;
        bool did_work = true;
        for(HalfEdgeID h: hevec)
            if(precond_collapse_edge(m, h)) {
                garbage.push_back(m.walker(h).opp().vertex());
                m.collapse_edge(h,false);
            }
            else {
                did_work = false;
                break;
            }
        
        if(did_work) {
            for(auto x: garbage)
                vset.erase(x);
            return true;
        }
        else m = m_safe;
    }
    return false;
}


bool remove_val2_vertices(Manifold& m, VertexSet& vset)
{
    int did_work = 0;
    vector<VertexID> garbage;
    for(VertexID v: vset)
        if(valency(m,v)==2)
        {
            Walker w = m.walker(v);
            if(precond_collapse_edge(m, w.halfedge())) {
                garbage.push_back(v);
                HalfEdgeID ph = w.prev().halfedge();
                FaceID f = w.prev().face();
                m.collapse_edge(w.halfedge());
                m.merge_faces(f, ph);
                ++ did_work;
            }
        }
    for(auto v: garbage)
        vset.erase(v);
    return did_work;
}


void console_spin_edge(MeshEditor* me, const std::vector<std::string> & args)
{
    Manifold& m = me->active_mesh();
    if(m.no_faces()< 10)
        return;
    me->save_active_mesh();
    auto sel_set = me->get_vertex_selection();
    HalfEdgeID spinner = InvalidHalfEdgeID;
    for(HalfEdgeID h: m.halfedges())
        if(m.in_use(h))
        {
            Walker w = m.walker(h);
            if(sel_set[w.vertex()] && sel_set[w.opp().vertex()])
            {
                spinner = h;
                break;
            }
        }
    if(spinner != InvalidHalfEdgeID) {
        Walker w = m.walker(spinner);
        sel_set[w.vertex()] = 0;
        sel_set[w.opp().vertex()] = 0;
        
        VertexID va = w.next().vertex();
        VertexID vb = w.opp().next().vertex();
        sel_set[va] = 1;
        sel_set[vb] = 1;
        me->get_vertex_selection() = sel_set;
        
        FaceID f = w.face();
        m.merge_faces(f, spinner);
        m.split_face_by_edge(f, va, vb);
    }
}

bool remove_gash(Manifold& m, FaceSet& fset)
{
    int did_work = 0;
    vector<FaceID> garbage;
    for(FaceID f: fset)
        if(m.in_use(f) && no_edges(m, f)==4)
        {
            array<VertexID,4> v3v = {InvalidVertexID, InvalidVertexID, InvalidVertexID, InvalidVertexID};
            array<int,4> val = {0,0,0,0};
            int k = 0;
            circulate_face_ccw(m, f, [&](VertexID v){
                val[k] = valency(m, v);
                v3v[k] = v;
                ++k;
            });
            int e_old = 0;
            for(int i=0;i<4;++i) e_old += sqr(val[i]-4);
            if(e_old>0) {
                int e[2] = {
                    (sqr(val[0]-1-4) +  sqr(val[2]-1-4) + sqr(val[1]+val[3]-2-4)),
                    (sqr(val[1]-1-4) +  sqr(val[3]-1-4) + sqr(val[0]+val[2]-2-4))
                };
                
                k= -1;
                if(e[0]<e_old && val[0]>=3 && val[2]>=3)
                    k = 1;
                if(e[1]<e[0] && val[1]>=3 && val[3]>=3)
                    k = 0;
                
                if(k>0)
                {
                    garbage.push_back(f);
                    FaceID fnew = m.split_face_by_edge(f, v3v[k], v3v[(k+2)%4]);
                    if(fnew != InvalidFaceID) {
                        Walker w = m.walker(fnew);
                        if(precond_collapse_edge(m, w.halfedge()))
                            m.collapse_edge(w.halfedge(), true);
                        ++did_work;
                    }
                    break;
                }
            }
        }
    for(auto f: garbage)
        fset.erase(f);
    return did_work;
}


void console_find_face_loops(MeshEditor* me, const std::vector<std::string> & args)
{
   
    Manifold& m = me->active_mesh();
    
    HalfEdgeAttributeVector<int> touched(m.no_vertices(), 0);
    for(auto f: m.faces()) {
        DebugRenderer::face_colors[f] = Vec3f(0);
    }
    for(auto h: m.halfedges())
        DebugRenderer::edge_colors[h] = Vec3f(0.4);
    for(auto v: m.vertices())
        DebugRenderer::vertex_colors[v] = Vec3f(1);

    int loop_no = 1;
    for(auto h: m.halfedges())
    {
        Walker w = m.walker(h);
        if(!touched[h]){
            vector<HalfEdgeID> hvec;
            set<HalfEdgeID> boundary;
            FaceAttributeVector<int> path_faces_touched(0);
            do {
                touched[w.halfedge()] = 1;
                touched[w.opp().halfedge()] = 1;
                if(path_faces_touched[w.face()] == 1)
                    DebugRenderer::face_colors[w.face()] = get_color(loop_no);
                DebugRenderer::edge_colors[w.halfedge()]= get_color(loop_no);
                DebugRenderer::edge_colors[w.opp().halfedge()]= get_color(loop_no);
                path_faces_touched[w.face()] = 1;
                hvec.push_back(w.halfedge());
                w = w.next().next().opp();
            }
            while(w.halfedge() != h && !touched[w.halfedge()]);
            ++ loop_no;
        }
    }

    
}

void console_quad_collapse(MeshEditor* me, const std::vector<std::string> & args)
{
    Manifold& m = me->active_mesh();
    me->save_active_mesh();
    auto sel_set = me->get_vertex_selection();
    for(HalfEdgeID h: m.halfedges())
        if(m.in_use(h))
        {
            Walker w = m.walker(h);
            if(sel_set[w.vertex()] && sel_set[w.next().next().vertex()])
            {
                FaceID f = w.face();
                VertexID v0 = w.vertex();
                VertexID v1 = w.next().next().vertex();
                FaceID fnew = m.split_face_by_edge(f, v0, v1);
                if(fnew != InvalidFaceID) {
                    Walker w = m.walker(fnew);
                    if(precond_collapse_edge(m, w.halfedge())) {
                        m.collapse_edge(w.halfedge(), true);
                        break;
                    }
                }
            }
        }
    console_find_face_loops(me, args);
}

void console_remove_valence2(MeshEditor* me, const std::vector<std::string> & args)
{
    Manifold& m = me->active_mesh();
    me->save_active_mesh();
    VertexSet verts;
    for(auto v: m.vertices()) verts.insert(v);
    while(remove_val2_vertices(m, verts));
}


void console_kill_loop(MeshEditor* me, const std::vector<std::string> & args)
{
    int doquit=0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >>  doquit;
    }

    Manifold& m = me->active_mesh();
    if(m.no_faces()<=6)
        return;
    
    me->save_active_mesh();

//    
//    VertexSet verts;
//    for(auto v: m.vertices()) verts.insert(v);
//    while(remove_val2_vertices(m, verts));

    
    HalfEdgeAttributeVector<int> touched(m.no_vertices(), 0);
    set<FaceID> fset_min;
    double fs_area_min = 1e200;
    set<HalfEdgeID> boundary_min;
    vector<HalfEdgeID> hvec_min;
    double asum_min = DBL_MAX;
    
    FaceAttributeVector<float> face_areas;
    for(auto f: m.faces()) {
     face_areas[f] = area(m,f);
        DebugRenderer::face_colors[f] = Vec3f(0);
    }
    for(auto h: m.halfedges())
        DebugRenderer::edge_colors[h] = Vec3f(0.4);
    for(auto v: m.vertices())
        DebugRenderer::vertex_colors[v] = Vec3f(1);
    
    for(auto h: m.halfedges())
    {
        Walker w = m.walker(h);
        if(!touched[h]){
            double asum = 0;
            vector<HalfEdgeID> hvec;
            set<HalfEdgeID> boundary;
            FaceAttributeVector<int> path_faces_touched(0);
            do {
                if(valency(m, w.vertex())+valency(m, w.opp().vertex())==6) break;
                touched[w.halfedge()] = 1;
                path_faces_touched[w.face()] = 1;
                asum += length(m, w.halfedge());
                hvec.push_back(w.halfedge());
                auto bh = w.prev().opp().halfedge();
                boundary.insert(bh);
                w = w.next().next().opp();
            }
            while(w.halfedge() != h && !touched[w.halfedge()] && !path_faces_touched[w.face()]);
            
            FaceAttributeVector<int> f_touched(m.allocated_faces(), 0);
            if(w.halfedge() == h) {
                double fs_area = 0;
                set<FaceID> face_set;
                queue<FaceID> fq;
                for(HalfEdgeID h: boundary) {
                    Walker w_bdry = m.walker(h);
                    if(boundary.count(w_bdry.opp().halfedge())==0) {
                        FaceID f_bdry = w_bdry.face();
                        if(!f_touched[f_bdry]) {
                            fq.push(f_bdry);
                            face_set.insert(f_bdry);
                            f_touched[f_bdry] = 1;
                            fs_area += face_areas[f_bdry];
                        }
                    }
                }
                while(!fq.empty())
                {
                    FaceID f = fq.front();
                    fq.pop();
                    circulate_face_ccw(m, f, [&](Walker w){
                        if(!boundary.count(w.halfedge())) {
                            FaceID fo = w.opp().face();
                            if(!f_touched[fo]) {
                                fq.push(fo);
                                face_set.insert(fo);
                                fs_area += face_areas[fo];
                                f_touched[fo] = 1;
                            }
                        }
                    });
                }
                if(fs_area< fs_area_min && hvec.size()>0) {
                    boundary_min = boundary;
                    fs_area_min = fs_area;
                    fset_min = face_set;
                    asum_min = asum;
                    hvec_min = hvec;
                }
            }
        }
    }
    if(boundary_min.size()>0) {
        //cout << "killing " << fs_area_min << " " << fss_min << endl;
        
        if(doquit) {
            for(auto f: fset_min)
                DebugRenderer::face_colors[f] = Vec3f(1,0,0);
            
            for(auto h: boundary_min)
            {
                Walker w = m.walker(h);
                DebugRenderer::edge_colors[h] = Vec3f(1,0,1);
                DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(1,0,1);
                
            }
            
            for(auto h: hvec_min)
            {
                Walker w = m.walker(h);
                DebugRenderer::edge_colors[h] = Vec3f(0,1,0);
                DebugRenderer::edge_colors[w.opp().halfedge()] = Vec3f(0,1,0);
                
            }
            return;
        }

        for(int iter=0;iter<3;++iter)
        for(auto h: hvec_min) if (m.in_use(h))
        {
            if(precond_collapse_edge(m, h))
                m.collapse_edge(h,true);
//            else {
////                me->restore_active_mesh();
//                for(auto h: m.halfedges())
//                    DebugRenderer::edge_colors[h] = Vec3f(0);
//                DebugRenderer::edge_colors[h] = Vec3f(1);
//                return;
//            }
        }
        
        for(auto h: hvec_min) if (m.in_use(h))
        {
            me->restore_active_mesh();
            return;
        }
        
        //while(remove_gash(m, fset_min));
        VertexSet vertex_set;
        for(auto f: fset_min)
            circulate_face_ccw(m, f, [&](VertexID v){vertex_set.insert(v);});
        for(auto h:boundary_min) {
            Walker w = m.walker(h);
            vertex_set.erase(w.vertex());
        }
//
//        while(remove_wedge(m, vertex_set));
//        while(remove_val2_vertices(m, vertex_set));
        
        for(int iter=0;iter<10;++iter)
        {
            for(auto v: vertex_set)
            {
                Vec3d avg_pos(0);
                avg_pos /= circulate_vertex_ccw(m,v, [&](VertexID vn){avg_pos += m.pos(vn);});
                m.pos(v) = avg_pos;
            }
        }
        
       if(!doquit)
           m.cleanup();
    }
    
}


void console_smooth_geodesic(MeshEditor* me, const std::vector<std::string> & args)
{
    int m1=1;
    int m2=2;
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> m1;
    }
    if(args.size() > 1){
        istringstream a0(args[1]);
        a0 >> m2;
    }
    if(args.size() > 2){
        istringstream a0(args[2]);
        a0 >> iter;
    }
    
    
    double alpha = 0.5;
    if(args.size() > 3){
        istringstream a0(args[3]);
        a0 >> alpha;
    }
    
    me->save_active_mesh();
    smooth_geodesic(me->get_mesh(m1-1) , me->get_mesh(m2-1), iter, alpha);
}

void console_refit(MeshEditor* me, const std::vector<std::string> & args)
{
    int m1=1;
    int m2=2;
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> m1;
    }
    if(args.size() > 1){
        istringstream a0(args[1]);
        a0 >> m2;
    }
    if(args.size() > 2){
        istringstream a0(args[2]);
        a0 >> iter;
    }
    
    
    double alpha = 0.5;
    bool preserve_poles = true;
    if(args.size() > 3){
        istringstream a0(args[3]);
        a0 >> alpha;
    }
    
    if(args.size() > 4){
        istringstream a0(args[4]);
        a0 >> preserve_poles;
    }
    
    me->save_active_mesh();
    smooth_and_refit(me->get_mesh(m1-1) , me->get_mesh(m2-1), iter, alpha, preserve_poles);
}

void console_remesh(MeshEditor* me, const std::vector<std::string> & args)
{
    int m1=1;
    int m2=2;
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> m1;
    }
    if(args.size() > 1){
        istringstream a0(args[1]);
        a0 >> m2;
    }
    if(args.size() > 2){
        istringstream a0(args[2]);
        a0 >> iter;
    }
    
    me->save_active_mesh();
    for(int i=0;i<iter;i+=10) {
        simplify_polar_mesh(me->active_mesh(), 0, 10);
        smooth_and_refit(me->get_mesh(m1-1) , me->get_mesh(m2-1), 3, 0.5, true);
    }
}



VertexAttributeVector<double> harmonic_extrema;
VertexAttributeVector<double> harmonic;

void console_save_assigned_vertex_values(MeshEditor* me, const std::vector<std::string> & args) {
    string str = me->active_visobj().file_name() + "-assigned.txt";
    ofstream ofs(str.c_str());
    for(VertexID v: me->active_mesh().vertices())
        ofs << harmonic_extrema[v] << " " << harmonic[v] << endl;
}

void console_load_assigned_vertex_values(MeshEditor* me, const std::vector<std::string> & args) {
    harmonic_extrema = VertexAttributeVector<double>(me->active_mesh().allocated_vertices(),0);
    string str = me->active_visobj().file_name() + "-assigned.txt";
    ifstream ifs(str.c_str());
    for(VertexID v: me->active_mesh().vertices())
        ifs >> harmonic_extrema[v] >> harmonic[v];
    me->active_visobj().get_scalar_field_attrib_vector() = harmonic;
    
}


void console_assign_values_to_vertex(MeshEditor* me, const std::vector<std::string> & args)
{
    double val=0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> val;
    }
    
    for(VertexID v: me->active_mesh().vertices())
        if(me->get_vertex_selection()[v]==1)
            harmonic_extrema[v] = val;
}

void console_quadrangulate_pole(MeshEditor* me, const std::vector<std::string> & args)
{
    me->save_active_mesh();
    pole_quadrangulation(me->active_mesh());
}

void console_heat_kernel_curvature(MeshEditor* me, const std::vector<std::string> & args)
{
    double R=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> R;
    }
    
    me->active_visobj().get_scalar_field_attrib_vector() = heat_kernel_laplacian(me->active_mesh(), R);
}


//void console_edge_equalize(MeshEditor* me, const std::vector<std::string> & args)
//{
//    int m1=1;
//    int m2=2;
//    int iter=1;
//    if(args.size() > 0){
//        istringstream a0(args[0]);
//        a0 >> m1;
//    }
//    if(args.size() > 1){
//        istringstream a0(args[1]);
//        a0 >> m2;
//    }
//    if(args.size() > 2){
//        istringstream a0(args[2]);
//        a0 >> iter;
//    }
//
//
//    double alpha = 0.5;
//    if(args.size() > 3){
//        istringstream a0(args[3]);
//        a0 >> alpha;
//    }
//
//    Manifold& m = me->get_mesh(m1-1);
//    Manifold& m_ref = me->get_mesh(m2-1);
//
//
//
//    vector<float> edge_lengths;
//    int n=0;
//    for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h,++n)
//        edge_lengths.push_back(length(m, *h));
//    sort(edge_lengths.begin(), edge_lengths.end());
//    float target_length = alpha*edge_lengths[n/2];
//
//    smooth_geodesic(m, me->get_mesh(m2-1), 20, 0.5);
//    for(int _iter=0;_iter<iter;++_iter)
//    {
//        vector<HalfEdgeID> short_edges, long_edges;
//        for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h)
//            if(*h<m.walker(*h).opp().halfedge())
//            {
//                float l = length(m,*h);
//                if(l> (4/3.0) * target_length)
//                    long_edges.push_back(*h);
//                else if(l< (4/5.0) * target_length)
//                    short_edges.push_back(*h);
//            }
//        for(int i=0;i<long_edges.size(); ++i)
//            if(m.in_use(long_edges[i]))
//                m.split_edge(long_edges[i]);
//        shortest_edge_triangulate(m);
//
//        for(int i=0;i<short_edges.size(); ++i)
//            if(m.in_use(short_edges[i]) && precond_collapse_edge(m, short_edges[i]))
//            {
//                Walker w = m.walker(short_edges[i]);
//                Vec3d mid_point = 0.5*(m.pos(w.vertex())+m.pos(w.opp().vertex()));
//                bool illegal = false;
//                for(Walker wc=m.walker(w.vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
//                    if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* target_length)
//                        illegal = true;
//                for(Walker wc=m.walker(w.opp().vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
//                    if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* target_length)
//                        illegal = true;
//
//                if(!illegal)
//                {
//                    if(boundary(m,w.vertex()) && !boundary(m,w.opp().vertex()))
//                        m.collapse_edge(short_edges[i],false);
//                    else if(!boundary(m,w.vertex()) && boundary(m,w.opp().vertex()));
//                    else
//                        m.collapse_edge(short_edges[i],true);
//                }
//            }
//        // optimize_valency(m);
//        //        minimize_dihedral_angle(m);
//
//        m.cleanup();
//
//
//
//
//        smooth_geodesic(m, me->get_mesh(m2-1), 20, 0.5);
//
//        maximize_min_angle(m, 0.97);
//    }
//}

#define GEODESIC_SPLIT_COLLAPSE

void console_edge_equalize(MeshEditor* me, const std::vector<std::string> & args)
{
    int m1=1;
    int m2=2;
    int outer_iter=5;
    int smooth_iter=10;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> m1;
    }
    if(args.size() > 1){
        istringstream a0(args[1]);
        a0 >> m2;
    }
    if(args.size() > 2){
        istringstream a0(args[2]);
        a0 >> outer_iter;
    }
    if(args.size() > 3){
        istringstream a0(args[3]);
        a0 >> smooth_iter;
    }
    double alpha = 0.5;
    if(args.size() > 4){
        istringstream a0(args[4]);
        a0 >> alpha;
    }
    
    Manifold& m = me->get_mesh(m1-1);
    Manifold& m_ref = me->get_mesh(m2-1);
    
    vector<float> edge_lengths;
    int n=0;
    for(HalfEdgeIDIterator h = m.halfedges_begin(); h != m.halfedges_end();++h,++n)
        edge_lengths.push_back(length(m, *h));
    sort(edge_lengths.begin(), edge_lengths.end());
    float median_length = edge_lengths[n/2];
    float _target_length = alpha*median_length;
    
    smooth_geodesic(m, me->get_mesh(m2-1), smooth_iter, 0.5);
    for(int _iter=0;_iter<outer_iter;++_iter)
    {
        float interp = _iter/float(outer_iter-1);
        float target_length = interp * _target_length + (1-interp) * median_length;
        cout << "Iter " << _iter << " target len " << target_length << endl;
#ifdef GEODESIC_SPLIT_COLLAPSE
        vector<pair<HalfEdgeID, ManiPoint>> short_edges, long_edges;
        VertexAttributeVector<ManiPoint> ptav = register_faces(m, m_ref);
#else
        vector<pair<HalfEdgeID, Vec3d>> short_edges, long_edges;
#endif
        
        for(HalfEdgeID h : m.halfedges())
            if(h<m.walker(h).opp().halfedge())
            {
                float len_e3 = length(m,h);
                Walker w = m.walker(h);
#ifdef GEODESIC_SPLIT_COLLAPSE
                LogMap log_map(m_ref, ptav[w.vertex()], 5 * len_e3);
                
                Vec2f uv = log_map.barycentric_uv(ptav[w.vertex()]);
                Vec2f uvo =log_map.barycentric_uv(ptav[w.opp().vertex()]);
                
                float l = length(uv-uvo);
                ManiPoint mid_point = log_map.find_uv(0.5*(uv+uvo));
#else
                float l = len_e3;
                Vec3d mid_point = m.pos(w.vertex());
#endif
                if(l> (4/3.0) * target_length)
                    long_edges.push_back(make_pair(h,mid_point));
                
                else if(l< (4/5.0) * target_length)
                    short_edges.push_back(make_pair(h,mid_point));
                
            }
        cout << "short " << short_edges.size() << " long " << long_edges.size() << endl;
        random_shuffle(begin(short_edges), end(short_edges));
        random_shuffle(begin(long_edges), end(long_edges));
        
        cout << "Splitting long edges" << endl;
        
        if(long_edges.size()>0) {
            for(int i=0;i<long_edges.size(); ++i)
                if(m.in_use(long_edges[i].first)) {
                    VertexID v = m.split_edge(long_edges[i].first);
#ifdef GEODESIC_SPLIT_COLLAPSE
                    ptav[v] = long_edges[i].second;
#endif
                }
#ifdef GEODESIC_SPLIT_COLLAPSE
            manipoints_to_3D(m, m_ref, ptav);
#endif
            shortest_edge_triangulate(m);
        }
        cout << "Collapsing short edges" << endl;
        if(short_edges.size()>0) {
            for(int i=0;i<short_edges.size(); ++i)
                if(m.in_use(short_edges[i].first) && precond_collapse_edge(m, short_edges[i].first))
                {
                    Walker w = m.walker(short_edges[i].first);
#ifdef GEODESIC_SPLIT_COLLAPSE
                    LogMap log_map(m_ref, ptav[w.vertex()], 3 * target_length);
                    Vec2f mid_point = log_map.barycentric_uv(short_edges[i].second);
#else
                    Vec3d mid_point = short_edges[i].second;
#endif
                    bool illegal = false;
                    for(Walker wc=m.walker(w.vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
#ifdef GEODESIC_SPLIT_COLLAPSE
                        if(length(log_map.uv_coords(wc.vertex())-mid_point) > (4/3.0)* target_length)
#else
                            if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* target_length)
#endif
                                illegal = true;
                    for(Walker wc=m.walker(w.opp().vertex()); !wc.full_circle(); wc = wc.circulate_vertex_ccw())
#ifdef GEODESIC_SPLIT_COLLAPSE
                        if(length(log_map.uv_coords(wc.vertex())-mid_point) > (4/3.0)* target_length)
#else
                            if(length(m.pos(wc.vertex())-mid_point) > (4/3.0)* target_length)
#endif
                                illegal = true;
                    
                    if(!illegal) {
#ifdef GEODESIC_SPLIT_COLLAPSE
                        VertexID v = w.vertex();
                        m.collapse_edge(w.halfedge(),true);
                        ptav[v] = short_edges[i].second;
#else
                        m.collapse_edge(w.halfedge(),true);
#endif
                    }
                }
#ifdef GEODESIC_SPLIT_COLLAPSE
            manipoints_to_3D(m, m_ref, ptav);
#endif
        }
        cout << "Optimizing" << endl;
        optimize_valency(m);
        cout << "Smoothing" << endl;
        smooth_geodesic(m, me->get_mesh(m2-1), smooth_iter, 0.5);
        cout << "Done" << endl;
    }
    m.cleanup();
}


void console_peek_values_at_vertex(MeshEditor* me, const std::vector<std::string> & args)
{
    double seleceted_extrema = 0.0;
    for(VertexID v: me->active_mesh().vertices())
        if(me->get_vertex_selection()[v]==1) {
            me->printf("harmonic, set value : %f, %f", harmonic[v], seleceted_extrema = harmonic_extrema[v]);
            me->printf("id %f", DebugRenderer::face_colors[me->active_mesh().walker(v).face()][2]);
        }
    if(seleceted_extrema != 0.0)
    {
        for(VertexID v: me->active_mesh().vertices())
            if(harmonic_extrema[v] == seleceted_extrema)
                me->get_vertex_selection()[v]=1;
        
    }
}

void console_reset_vertex_values(MeshEditor* me, const std::vector<std::string> & args)
{
    harmonic_extrema = VertexAttributeVector<double>(me->active_mesh().allocated_vertices(),0);
}


void console_skeletonize(MeshEditor* me, const std::vector<std::string> & args)
{
    double t=0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
    }
    me->save_active_mesh();
    //    skeletonize(me->active_mesh(),me->active_visobj().get_scalar_field_attrib_vector(), t);
    skeleton_line_field(me->active_mesh(), me->active_visobj().get_scalar_field_attrib_vector(), me->active_visobj().get_line_field_attrib_vector());
}

void console_make_harmonic(MeshEditor* me, const std::vector<std::string> & args)
{
    double t=0;
    bool doshift = false;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
        doshift = true;
    }
    me->save_active_mesh();
    
    VertexAttributeVector<int> nailed(me->active_mesh().allocated_vertices(),0);
    harmonic.resize(me->active_mesh().allocated_vertices());
    for(VertexID v: me->active_mesh().vertices()) {
        if(harmonic_extrema[v]!= 0.0)
            nailed[v] = 1;
        harmonic[v] = harmonic_extrema[v];
    }
    smooth_fun(me->active_mesh(), nailed,
               harmonic, 10000);
    
    if (doshift) {
        for(VertexID v: me->active_mesh().vertices())
            harmonic[v] += t;
        
    }
    
    me->active_visobj().get_scalar_field_attrib_vector() = harmonic;
}

void console_constrained_smooth(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1000;
    bool doshift = false;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    auto P = me->active_mesh().positions_attribute_vector();
    smooth_fun(me->active_mesh(), me->get_vertex_selection(),
               P, iter);
    me->active_mesh().positions_attribute_vector() = P;
}


void console_make_adf(MeshEditor* me, const std::vector<std::string> & args)
{
    double t=0;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
    }
    double vmin, vmax;
    make_adf_fun(me->active_mesh(), t, harmonic, vmin, vmax);
    me->active_visobj().get_scalar_field_attrib_vector() = harmonic;
    
}
void console_make_height_fun(MeshEditor* me, const std::vector<std::string> & args)
{
    double t=0;
    bool doshift = false;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> t;
        doshift = true;
    }
    double vmin, vmax;
    make_height_fun(me->active_mesh(), harmonic, vmin, vmax);
    if (doshift) {
        for(VertexID v: me->active_mesh().vertices())
            harmonic[v] += t;
        
    }
    me->active_visobj().get_scalar_field_attrib_vector() = harmonic;
    
}


void console_laplacian_flatten(MeshEditor* me, const std::vector<std::string> & args)
{
    int iter=1;
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> iter;
    }
    me->save_active_mesh();
    Manifold& m = me->active_mesh();
    auto laplacian = [&](VertexID v)
    {
        Vec3d p(0);
        int n = circulate_vertex_ccw(m, v, (std::function<void(VertexID)>)[&](VertexID v){ p += m.pos(v); });
        return p / n - m.pos(v);
    };
    VertexAttributeVector<Vec3d> new_pos = m.positions_attribute_vector();
    static VertexAttributeVector<Vec3d> pdir;
    static bool was_here = false;
    if(!was_here) {
        was_here = true;
        skeleton_line_field(m, me->active_visobj().get_scalar_field_attrib_vector(), pdir);
    }
    VertexAttributeVector<Vec3d> safe_pos = m.positions_attribute_vector();
    
    //    VertexAttributeVector<Vec3d> normal_field(m.allocated_vertices(),Vec3d(0));
    //    for(VertexID v: m.vertices())
    //        normal_field[v] = normal(m, v);
    //    Harmonics h(m);
    //    auto normal_field_t = h.analyze_signal(normal_field);
    //    pdir = h.reconstruct_signal(normal_field_t,2);
    
    for(int i=0;i<iter;++i) {
        for(VertexID v: m.vertices())
        {
            Vec3d dir = pdir[v];//normalize(global_dir - pdir[v]*dot(global_dir, pdir[v]));
            new_pos[v] = m.pos(v) + 0.5 * dir*dot(dir,laplacian(v));
            Vec3d n = normal(m, v);
            
        }
        m.positions_attribute_vector() = new_pos;
    }
    for(VertexID v: m.vertices())
        me->active_visobj().get_scalar_field_attrib_vector()[v] = dot(normal(m,v),pdir[v]);
    HMesh::VertexAttributeVector<int> nailed(m.no_vertices(),0);
    smooth_fun(m, nailed, me->active_visobj().get_scalar_field_attrib_vector(), 1);
    //  m.positions_attribute_vector() = safe_pos;
}





void load_labels(MeshEditor* me, const std::vector<std::string> & args)
{
    Manifold& m =  me->active_mesh();
    for(auto f: m.faces())
        DebugRenderer::face_colors[f] = Vec3f(0);
    
    string label_path = "/Users/jab/Studio/3DModels/SegmentationBenchmark/results/";
    const string& path = me->file_name();
    
    smatch match;
    regex_search(path, match, regex("(\\w+)\\.off"));
    
    string seg_file = label_path + match[1].str() + ".lab";
    
    ifstream f(seg_file.data());
    
    vector<FaceID> faces;
    for(auto fid: me->active_mesh().faces())
        faces.push_back(fid);
    
    int n = 0;
    while(!f.eof() && f.good())
    {
        int i;
        
        while (f>>i) {
            if(i-1<faces.size()) {
                DebugRenderer::face_colors[faces[i-1]] = get_color(n);
            }
        }
        if(!f.eof())
            f.clear();
        string str;
        f >> str;
        cout << str << endl;
        ++n;
        
    }
    
    for(auto h: m.halfedges())
    {
        Walker w = m.walker(h);
        if(DebugRenderer::face_colors[w.face()]!=DebugRenderer::face_colors[w.opp().face()])
            DebugRenderer::edge_colors[h] = Vec3f(0);
        else
            DebugRenderer::edge_colors[h] = 0.5 * DebugRenderer::face_colors[w.face()];
    }
    for(auto v: m.vertices())
        DebugRenderer::vertex_colors[v] = DebugRenderer::edge_colors[m.walker(v).halfedge()];
    
}

void console_polarize(MeshEditor* me, const std::vector<std::string> & args)
{
    int divisions = 50;
    
    if(args.size() > 0){
        istringstream a0(args[0]);
        a0 >> divisions;
    }
    
    double t=1;
    
    if(args.size() > 1){
        istringstream a0(args[1]);
        a0 >> t;
    }
    
    me->save_active_mesh();
    
    double vmin=DBL_MAX, vmax=-DBL_MAX;
    for(VertexID v: me->active_mesh().vertices()) {
        vmin = min(harmonic[v], vmin);
        vmax = max(harmonic[v], vmax);
    }
    polarize_mesh(me->active_mesh(), harmonic, vmin, vmax, divisions);
    
}


void register_console_funcs(MeshEditor* me)
{
    me->register_console_function("polar.pole_quadrangulation", console_quadrangulate_pole,
                                  "polar.pole_quadrangulation");
    
    me->register_console_function("polar.load_assigned_vertex_values", console_load_assigned_vertex_values,
                                  "polar.load_assigned_vertex_values");
    me->register_console_function("polar.save_assigned_vertex_values", console_save_assigned_vertex_values,
                                  "polar.save_assigned_vertex_values");
    me->register_console_function("polar.peek", console_peek_values_at_vertex,
                                  "polar.peek");
    me->register_console_function("polar.reset_vertex_values", console_reset_vertex_values,
                                  "polar.reset_vertex_values");
    me->register_console_function("polar.make_adf", console_make_adf,
                                  "polar.make_adf");
    me->register_console_function("polar.make_height_fun", console_make_height_fun,
                                  "polar.make_height_fun");
    
    me->register_console_function("polar.make_harmonic", console_make_harmonic,
                                  "polar.make_harmonic");
    me->register_console_function("polar.assign_vertex_values", console_assign_values_to_vertex,
                                  "polar.assign_vertex_values");
    
    me->register_console_function("polar.convert", console_polarize,
                                  "polar.convert <cuts> <t>");
    me->register_console_function("polar.segment", console_polar_segment,
                                  "polar.segments <show face segments>");
    me->register_console_function("polar.skeleton", console_polar_skeleton,
                                  "polar.skeleton <fraction>");
    me->register_console_function("polar.show_skin", console_polar_show_skin,
                                  "polar.show_skin <joint>");
    
    me->register_console_function("polar.subdivide", console_polar_subdivide,
                                  "polar.subdivide <iterations>");
    me->register_console_function("polar.simplify", console_simplify_polar,
                                  "polar.simplify <frac> <iter>");
    me->register_console_function("polar.branch", console_polar_add_branch,
                                  "polar.branch");
    me->register_console_function("polar.refine_poles", console_refine_poles,
                                  "polar.refine_poles");
    me->register_console_function("polar.refit", console_refit,
                                  "polar.refit");
    me->register_console_function("polar.remesh", console_remesh,
                                  "polar.remesh");
    me->register_console_function("polar.dist_color", console_dist_color,
                                  "polar.dist_color");
    me->register_console_function("smooth.geodesic", console_smooth_geodesic,
                                  "smooth.geodesic");
    me->register_console_function("optimize.edge_equalize", console_edge_equalize,
                                  "optimize.edge_equalize");
    me->register_console_function("heat_kernel_curvature", console_heat_kernel_curvature,
                                  "heat_kernel_curvature");
    me->register_console_function("load_segmentation",
                                  load_labels,
                                  "");
    
    me->register_console_function("spin_edge", console_spin_edge, "");
    me->register_console_function("find_loops", console_find_face_loops, "");
    me->register_console_function("kill_loop", console_kill_loop, "");
    me->register_console_function("kill_v2", console_remove_valence2, "");
    me->register_console_function("quad_collapse", console_quad_collapse, "");
    me->register_console_function("polar.skeletonize", console_skeletonize, "");
    me->register_console_function("polar.smooth_constrained", console_constrained_smooth, "");
    
    
    me->register_console_function("planarize.flatten", console_laplacian_flatten, "");
    me->register_console_function("polar.collapse_selected", console_polar_collapse_selected, "");
    me->register_console_function("polar.refine_selected", console_polar_refine_selected, "");
    me->register_console_function("polar.merge_poles", console_merge_poles, "");
}
