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

#include <map>
#include <queue>
#include <algorithm>

#include "Procedural/Helpers/Plane.h"
#include <MeshEditE/Procedural/Operations/geometric_operations.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>

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


