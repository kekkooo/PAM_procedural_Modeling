//
//  Test.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 25/07/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "Test.h"
#include <GEL/GLGraphics/MeshEditor.h>
#include <GEL/CGLA/Vec3d.h>
#include <GEL/HMesh/obj_save.h>
#include <GEL/CGLA/Mat3x3d.h>

#include <map>
#include <queue>
#include <algorithm>

#include "Procedural/Helpers/Plane.h"
#include <MeshEditE/Procedural/Operations/geometric_operations.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <MesheditE/Procedural/Helpers/structural_helpers.h>
#include "patch_mapping.h"


using namespace GLGraphics;
using namespace std;
using namespace CGLA;
using namespace HMesh;
using namespace Procedural::Geometry;
using namespace Procedural::Operations::Geometric;


bool test_ring_barycenter( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    vector< VertexID > vs;
    Vec3d b = ring_barycenter(m, h, vs);
    assert( vs.size() > 0 );
    Plane p( m.pos( vs.front()), m.pos(vs.back()), m.pos(vs[ vs.size()/2 ]));
    
    for( VertexID v : vs )
    {
        assert( p.OnPlane( m.pos( v )));
    }
    assert( p.OnPlane( b));

    return true;
}


Vec3d simple_random_direction ( Manifold& m, VertexID v, double vertex_normal_weight )
{
    // take the vertex normal
    // choose randomly the normal of one of the sorrounding faces
    // sum them
    // multiply for a random ration between the range [ -0.5, 0.5 ]
    Vec3d v_normal = vertex_normal( m, v );
    int     times   = rand() % 23;
    double  ratio   = ( (double)( rand() % 1000 ) - 500.0 )/10000.0;
//    cout << "times is : " << times << " ratio is : " << ratio <<endl;
    Walker w = m.walker(v);
    for( int i = 0; i < times; w = w.circulate_vertex_ccw(), i++ );
    
    Vec3d f_normal  = face_normal( m, w.face( ));
    
    f_normal.normalize();
    v_normal.normalize();
    
    Vec3d dir       = ( v_normal + f_normal ) * ratio;
    return dir;
}

Vec3d alt_simple_random_direction ( Manifold& m, VertexID   v )
{
    Vec3d v_normal = vertex_normal( m, v );
    int     times   = rand() % 23;
    double  ratio   = ( (double)( rand() % 1000 ) - 500.0 )/10000.0;
//    cout << "times is : " << times << " ratio is : " << ratio <<endl;
    Walker w = m.walker(v);
    for( int i = 0; i < times; w = w.circulate_vertex_ccw(), i++ );
    
    Vec3d other_dir  = m.pos( w.vertex()) - m.pos(v);
    
    other_dir.normalize();
    v_normal.normalize();
    
    Vec3d dir       = ( v_normal + other_dir * 1.5 ) * ratio;
    return dir;
}

Vec3f color_ramp( int value, int max )
{
    // la colorazione va da :
    // rosso    : molto rosso, poco verde,  poco blu
    // arancio  : molto rosso, verde medio, poco blu
    // verde    : poco  rosso, molto verde, poco blu

    assert( value >= 0  );
    assert( value <= max );
    assert( max > 1     );
    
    float   red     = 1.0,
            blue    = 0.1,
            green   = 0.1;
    
    int red_scale   = value;
    int green_scale = max - value;
    
    green   = green_scale / (float)max;
    red     = red_scale  / (float)max;
    
    assert( green >= 0 && green <=1 );
    assert( red >= 0 && red <=1 );
    
//    cout << value << " / " << max << endl << "bring to :" << red << ", " << green << ", " << blue << endl;
    
    return Vec3f( red, green, blue );
}

Vec3f color_ramp2( int value, int max )
{
    float er = 240.0,
          eg = 50.0,
          eb = 50.0,
          sr = 5.0,
          sg = 180.0,
          sb = 60;
    float n;
    float r, g, b;
    n = (float)value / (float)max;
    r = (float)sr * ( 1.0f -n ) + (float)er * n;
    g = (float)sg * ( 1.0f -n ) + (float)eg * n;
    b = (float)sb * ( 1.0f -n ) + (float)eb * n;
    
    return Vec3f( r/256.0, g/256.0, b/256.0 );
}

void linspace ( double min, double max, int num, std::vector<double> &values)
{
    if ( num <= 0 ) return;
    
    double pace = ( max - min ) / (double)( num - 1 );
    values.push_back( min );
    for(int i=1; i< num; i++)
    {
        values.push_back(values[i-1] + pace);
    }
}

CGLA::Mat4x4d get_rotation_mat4d ( CGLA::Vec3d axis, double angle )
{
    CGLA::Mat4x4d t;
    CGLA::Vec3d ax( 1.0, 0.0, 0.0 ), ay( 0.0, 1.0, 0.0 ), az( 0.0, 0.0, 1.0 );
    double cosine = cos(angle), sine = sin(angle);
    axis.normalize();

    double x = axis[0], y = axis[1], z = axis[2];
    
    CGLA::Vec4d row1( ( x * x ) + (( 1.0 - x * x ) * ( cosine )),
                      (( 1.0 - ( cosine )) * x * y) - (( sine ) * z ),
                      (( 1.0 - ( cosine )) * x * z) + (( sine ) * y ),
                      0.0 );
    CGLA::Vec4d row2( (( 1.0 - ( cosine)) * y * x) + (( sine ) * z),
                      ( y * y ) + (( 1.0 - y * y ) * (cosine)),
                      (( 1.0 - ( cosine)) * y * z) - (( sine ) * x), 0.0 );
    
    CGLA::Vec4d row3( (( 1.0 - ( cosine)) * z * x) - (( sine ) * y),
                      (( 1.0 - ( cosine)) * z * y) + (( sine ) * x),
                      ( z * z ) + (( 1.0 - z * z ) * (cosine)), 0.0 );
    CGLA::Vec4d row4( 0.0, 0.0, 0.0, 1.0 );

    t = CGLA::Mat4x4d( row1, row2, row3, row4 );
    return t;
}

CGLA::Mat4x4d get_alignment_for_2_vectors( Vec3d v1, Vec3d v2, Vec3d centroid )
{
    
    assert( !(( v1[0] > GEO_EPS ) && ( v1[1] > GEO_EPS ) && ( v1[2] > GEO_EPS )));
    assert( !(( v2[0] > GEO_EPS ) && ( v2[1] > GEO_EPS ) && ( v2[2] > GEO_EPS )));
    
    v1.normalize();
    v2.normalize();
    Vec3d   rotation_axis;
    double  rotation_angle  = get_angle( v1, v2 );
    
    if( dot ( v1, v2 ) > 0.8 ){
        rotation_axis = CGLA::cross( v1, v2 );
    }
    else{
        rotation_axis = CGLA::cross( v1, -v2 );
        rotation_angle = -rotation_angle;
    }
    
    bool parallel = ( fabs(rotation_axis[0] + rotation_axis[1] + rotation_axis[2]) < GEO_EPS );
    bool opposite = ( (( v1[0] + v2[0] ) < GEO_EPS ) && (( v1[1] + v2[1] ) < GEO_EPS ) && (( v1[2] + v2[2] ) < GEO_EPS ));
    
    if( parallel && opposite ){
        // probably v1 and v2 are equal
        cout << v1 << " is equal to " << v2 << "?" << endl;
        cout << "rotation axis :" << endl << rotation_axis << endl;
        cout << "angle :" << endl << rotation_angle << endl;
        cout << "build reflection" << endl;
        Mat4x4d T;
        v1.normalize();
//        buildReflectionMatrix(v1, T);
        // instead of building a reflection matrix
        // build a transformation that brings v1 be aligned with Y-axis and centroid to be placed at origin
        // then build a rotation matrix of angle PI along X-axis
        Vec3d   y_axis  = Vec3d(0.0, 1.0, 0.0);
        Vec3d   _axis   = CGLA::cross( v1, y_axis );
        double  _angle  = get_angle( v1, y_axis );
        Mat4x4d rot_to_y, rot;
        // normal could be aligned to y_axis
        if( _axis[0] + _axis[1] + _axis[2] < GEO_EPS ){
            rot_to_y = identity_Mat4x4d();
        }
        else{
            rot_to_y = get_rotation_mat4d( _axis, _angle);
        }
        
        
        Mat4x4d rot_to_y_invert =  invert_ortho(rot_to_y);
        
        // What happens if V1 is the XAxis or the ZAxis??
        rot = rotation_Mat4x4d( XAXIS, M_PI );
        Mat4x4d tr_origin = translation_Mat4x4d( -centroid );
        Mat4x4d tr_back   = translation_Mat4x4d( centroid );
        T = rot_to_y_invert * tr_back * rot * tr_origin * rot_to_y;

        cout << T << endl;
        if( isnan(T[1][1])){
            cout << tr_back << endl << rot << endl << tr_origin << endl << rot_to_y << endl;
        }
        
        
        assert( !isnan( T[1][1] ));
        return T;
    }
    
    // if angle is near 0 but normals are opposite just return identity
    if( parallel ){
        return identity_Mat4x4d();
    }

    
    assert( !isnan(rotation_angle )); // fail!
    
    Mat4x4d rot             = get_rotation_mat4d( rotation_axis, rotation_angle);
    // find center of the module mesh
    Mat4x4d tr_origin = translation_Mat4x4d( -centroid );
    Mat4x4d tr_back   = translation_Mat4x4d( centroid );
    Mat4x4d t =  tr_back * rot * tr_origin;

//    cout << rot         << endl << endl;
//    cout << tr_origin   << endl << endl;
//    cout << tr_back     << endl << endl;
//    cout << t           << endl << endl;
    assert( !isnan( t[1][1] ));
    return t;
}

void alt_glue_poles(Manifold& mani, VertexID vid0, VertexID vid1)
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

void bezier ( CGLA::Vec3d p0, CGLA::Vec3d p1, CGLA::Vec3d p2, int n, vector< CGLA::Vec3d >& points )
{
    points.clear();
    points.push_back(p0);
    // get the pace on the interval [0, 1]
    double pace = 1.0 / (double)( n - 1 );
    for( int i = 1; i <= n; ++i )
    {
        double t = (double)i * pace;
        CGLA::Vec3d pi = (( 1.0 - t ) * ( 1.0 -t ) * p0 )
                        + ( 2 * t *( 1 - t ) * p1 )
                        + ( t * t ) * p2;
        points.push_back( pi );
    }
    points.push_back(p2);
}

void save_intermediate_result ( HMesh::Manifold &m, const std::string &folder_path, int step_number )
{
    stringstream oss;
    oss << folder_path << "step" << step_number << ".obj";
    HMesh::obj_save( oss.str(), m );
}

void save_colored_obj ( HMesh::Manifold& m, string &path )
{
    // assume that path does not includes the file extension
    string obj_filename = path + ".obj";
    string mtl_filename = path + ".mtl";
    // need to build the patch structure
    HMesh::FaceAttributeVector<int> face_patch;
//    polar_extract_patches( m, face_patch );
    build_patches(m, face_patch);
    int max_patch = numeric_limits<int>::min();
    
    // *** OBJ ***
    ofstream os( obj_filename.data( ));
    assert( !os.bad() );
    
    VertexAttributeVector<int> vmap;
    int k = 0;
    
    // mtllib <filename>.mtl
    os << "mtllib "<< mtl_filename << endl;

    // - elenco dei vertici v <x> <y> <z>
    for(VertexIDIterator v = m.vertices_begin(); v != m.vertices_end(); ++v){
        Vec3d p = m.pos(*v);
        os << "v "<< p[0] << " " << p[1] << " " << p[2] << "\n";
        vmap[*v] = k++;
    }
    
    // - elenco delle facce su due righe
    //      usemtl PATCH_X
    //      elenco vertici
    for(FaceIDIterator f = m.faces_begin(); f != m.faces_end(); ++f){
        vector<int> verts;
        for(Walker w = m.walker(*f); !w.full_circle(); w = w.circulate_face_ccw()){
            int idx = vmap[w.vertex()];
            assert(static_cast<size_t>(idx) < m.no_vertices());
            // move subscript range from 0..size-1 to 1..size according to OBJ standards
            verts.push_back(idx + 1);
        }
        os << "usemtl PATCH_" << face_patch[*f] << endl;;
        os << "f ";
        for(size_t i = 0; i < verts.size() ; ++i){
            os << verts[i];
            if (i+1==verts.size())
                break;
            os << " ";
        }
        os<<endl;
        
        // keep track of the number of patches
        if( face_patch[*f] > max_patch ) max_patch = face_patch[*f];
    }
    os.close();
    
    
    // *** MTL ***
    ofstream mtl_os( mtl_filename.data( ));
    assert( !mtl_os.bad() );
    
    // for each patch of the model there will be 2 rows on the MTL file
    for( int i = 0; i < max_patch; ++i )
    {
        mtl_os << "newmtl PATCH_" << i << endl;
        mtl_os << "Kd 0.535156 0.824219 0.894531" << endl;
    }
    // newmtl PATCH_X // where X is the id of the patch
    // Kd 0.894531 0.535156 0.535156
    
    mtl_os.close();

}

int get_starter_offset( Manifold &m1, VertexID p1, Manifold &m2, VertexID p2 )
{
    assert( valency( m1, p1 ) == valency( m2, p2 ));
    
    std::vector<HalfEdgeID> ring_p1, ring_p2;
    std::vector<double>     angles_p1, angles_p2, diffs;
    std::vector<Vec3d>      vecs_p1, vecs_p2;
    Vec3d                   p1_pos = m1.pos( p1 ),
                            p2_pos = m2.pos( p2 );
    Walker                  w1 = m1.walker( p1 ),
                            w2 = m2.walker( p2 );
    for( ; w1.full_circle(); w1 = w1.circulate_vertex_ccw() )
    {
        ring_p1.push_back( w1.halfedge( ));
        vecs_p1.push_back( m1.pos( w1.vertex()) - p1_pos );
    }
    for( ; w2.full_circle(); w2 = w2.circulate_vertex_ccw() )
    {
        ring_p2.push_back( w1.halfedge( ));
        vecs_p2.push_back( m2.pos( w2.vertex()) - p2_pos );
    }
    
    for( int i = 1; i < ring_p1.size(); ++i )
    {
        angles_p1.push_back( get_angle(vecs_p1[i-1], vecs_p1[i] ));
        angles_p2.push_back( get_angle(vecs_p2[i-1], vecs_p2[i] ));
    }
    angles_p1.push_back( get_angle( vecs_p1[vecs_p1.size() - 1], vecs_p1[0] ));
    angles_p2.push_back( get_angle( vecs_p2[vecs_p1.size() - 1], vecs_p2[0] ));
#pragma message "scritto di getto, da verificare"
    for( int i = 0; i < ring_p1.size(); ++i ){
        for( int j = 0; j < ring_p1.size(); ++j )
        {
            int i1 = j,
                i2 = ( i + j ) % ring_p1.size();
            double diff_sum = 0.0;
            for( int k = 0; k < ring_p1.size(); ++k )
            {
                diff_sum += fabs( angles_p1[i1] - angles_p2[i2] );
            }
            diffs.push_back( diff_sum );
        }
    }
    // find the index at the difference is min
    int min_index = 0;
    for( int i = 1; i < ring_p1.size(); ++i )
    {
        if( diffs[i] < diffs[min_index] ) { min_index = i; }
    }
    return min_index;
}

void bridge_pole_one_rings(Manifold& mani, VertexID vid0, VertexID vid1)
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

void buildReflectionMatrix ( CGLA::Vec3d& planeNormal, CGLA::Mat4x4d &T ){
    planeNormal.normalize();
    
    double  a = planeNormal[0],
            b = planeNormal[1],
            c = planeNormal[2];

    Vec4d r0 = Vec4d( 1.0 - 2.0 * a * a , -2.0 * a * b      ,-2.0 * a * c       ,0.0 );
    Vec4d r1 = Vec4d( -2.0 * a * b      , 1.0 - 2 * b * b   ,-2.0 * b * c       ,0.0 );
    Vec4d r2 = Vec4d( -2.0 * a * c      , -2.0 * b * c      ,1.0 - 2.0 * c * c  ,0.0 );
    Vec4d r3 = Vec4d( 0.0               ,0.0                ,0.0                ,1.0 );
    
    T = Mat4x4d( r0, r1, r2, r3 );
}







