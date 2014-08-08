//
//  geometric_properties.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "geometric_properties.h"

using namespace CGLA;
using namespace std;
using namespace HMesh;

namespace Procedural{
    namespace Geometry{


Vec3d ring_barycenter ( HMesh::Manifold& m, HMesh::HalfEdgeID h )
{
    vector< VertexID > vs;
    return ring_barycenter(m, h, vs);
}

Vec3d ring_barycenter ( Manifold& m, HalfEdgeID h, vector< VertexID > &vertices   )
{
    assert(h != InvalidHalfEdgeID);
    
    Walker w = m.walker( h );
    Vec3d barycenter(0);
    
    while ( !w.full_circle( ))
    {
        vertices.push_back( w.vertex( ));
        barycenter += m.pos( w.vertex() );
        w = w.next().opp().next();
    }
    
    barycenter /= vertices.size();
    return barycenter;
}

int valence ( Manifold& m, VertexID v )
{
    Walker w = m.walker( v );
    for( ; !w.full_circle(); w = w.circulate_vertex_ccw());
    return w.no_steps();
}

bool is_singularity ( Manifold& m, VertexID v )
{
    return (! ( valence(m, v) == 4));
}
        

// this function is though to be used only for PAMs
Vec3d face_normal ( Manifold& m, FaceID f)
{
    assert( m.in_use( f ));
    
    Walker w = m.walker(f);
    vector< VertexID > verts;
    vector< Vec3d > pos;
    
    for ( ; !w.full_circle(); w = w.circulate_face_ccw( ))
    {
        verts.push_back( w.vertex( ));
        pos.push_back( m.pos( w.vertex( )));
    }
    
    // check if it is a quad or a triangle
    bool is_triangle = w.no_steps() == 3;
    if ( is_triangle )
    {
        auto n = cross( pos[1] - pos[0], pos[2] - pos[0] );
        //        n.normalize();
        return n;
    }
    else
    {
        auto n1 = cross( pos[1] - pos[0], pos[2] - pos[0] );
        auto n2 = cross( pos[2] - pos[0], pos[3] - pos[0] );
        auto n = ( n1 + n2 )/2.0;
        //        n.normalize();
        return n;
    }
}

Vec3d vertex_normal ( Manifold& m, VertexID v)
{
    Walker w = m.walker(v);
    Vec3d vn(0);
    for ( ; !w.full_circle(); w = w.circulate_vertex_ccw( ))
    {
        vn += face_normal(m, w.face());
    }
    
    vn / w.no_steps();
    //    vn.normalize();
    return vn;
}


}}