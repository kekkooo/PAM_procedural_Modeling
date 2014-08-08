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
    cout << "times is : " << times << " ratio is : " << ratio <<endl;
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
    cout << "times is : " << times << " ratio is : " << ratio <<endl;
    Walker w = m.walker(v);
    for( int i = 0; i < times; w = w.circulate_vertex_ccw(), i++ );
    
    Vec3d other_dir  = m.pos( w.vertex()) - m.pos(v);
    
    other_dir.normalize();
    v_normal.normalize();
    
    Vec3d dir       = ( v_normal + other_dir * 1.5 ) * ratio;
    return dir;

}

