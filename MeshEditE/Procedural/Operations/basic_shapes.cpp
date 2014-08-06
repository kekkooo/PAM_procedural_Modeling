//
//  basic_shapes.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 06/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#include "basic_shapes.h"
#include <GEL/CGLA/Vec3d.h>

using namespace std;
using namespace CGLA;

void create_cube_vertices( double edge_length, vector<Vec3d> &vertices )
{
    float plus  = edge_length / 2.0;
    float minus = -plus;
    
    // initialize points
    vertices[0] = Vec3d( minus, plus, plus );
    vertices[1] = Vec3d( minus, plus, minus );
    vertices[2] = Vec3d( plus, plus, minus );
    vertices[3] = Vec3d( plus, plus, plus );
    vertices[4] = Vec3d( minus, minus, plus );
    vertices[5] = Vec3d( minus, minus, minus );
    vertices[6] = Vec3d( plus, minus, minus );
    vertices[7] = Vec3d( plus, minus, plus );
}


namespace Procedural
{
namespace Operations
{

    void create_cube_vertices( double edge_length, vector<Vec3f> &vertices )
    {
        float plus  = edge_length / 2.0;
        float minus = -plus;
        
        // initialize points
        vertices[0] = Vec3f( minus, plus, plus );
        vertices[1] = Vec3f( minus, plus, minus );
        vertices[2] = Vec3f( plus, plus, minus );
        vertices[3] = Vec3f( plus, plus, plus );
        vertices[4] = Vec3f( minus, minus, plus );
        vertices[5] = Vec3f( minus, minus, minus );
        vertices[6] = Vec3f( plus, minus, minus );
        vertices[7] = Vec3f( plus, minus, plus );
    }
    
    // Creates a cube with little piramids (with squared foot) instead of top and bottom face
    // One of the simplest PAMs
    void create_basic_PAM( HMesh::Manifold& m, double edge_length )
    {
        vector< Vec3f > vertices( 10 );
        vector< int >   faces   ( 12, 3 );
        vector< int >   indices ( 40 );
        
        create_cube_vertices(edge_length, vertices);
        // add the pole vertices;
        vertices[8] = Vec3f(0.0, ( 2.0 * edge_length / 3.0 ), 0.0 );
        vertices[9] = Vec3f(0.0, ( -2.0 * edge_length / 3.0 ), 0.0 );
        
        faces[0] = 4;   faces[1] = 4;   faces[2] = 4;   faces[3] = 4;
        
        // LEFT
        indices[0]  = 1;    indices[1]  = 5;    indices[2]  = 4;    indices[3]  = 0;
        // BACK
        indices[4]  = 1;    indices[5]  = 2;    indices[6]  = 6;    indices[7]  = 5;
        // FRONT
        indices[8]  = 3;    indices[9]  = 0;    indices[10] = 4;    indices[11] = 7;
        // RIGHT
        indices[12] = 7;    indices[13] = 6;    indices[14] = 2;    indices[15] = 3;
        // upper pole
        indices[16] = 2;    indices[17] = 1;    indices[18] = 8;
        indices[19] = 1;    indices[20] = 0;    indices[21] = 8;
        indices[22] = 0;    indices[23] = 3;    indices[24] = 8;
        indices[25] = 3;    indices[26] = 2;    indices[27] = 8;
        // lower pole
        indices[28] = 7;    indices[29] = 4;    indices[30] = 9;
        indices[31] = 4;    indices[32] = 5;    indices[33] = 9;
        indices[34] = 5;    indices[35] = 6;    indices[36] = 9;
        indices[37] = 6;    indices[38] = 7;    indices[39] = 9;
        
        m.build( vertices.size( ), reinterpret_cast< float * >( &vertices[0] ),
                faces.size( ), &faces[0], &indices[0] );
    }
    
    // Creates the triangle mesh of a cube
    void create_triangle_cube( HMesh::Manifold& m, double edge_length )
    {
        vector< Vec3f > vertices( 8 );
        vector< int >   faces   ( 12, 3 );
        vector< int >   indices ( 36 );
        
        create_cube_vertices(edge_length, vertices);
        
        // for each face save the indices of its vertices
        // UP
        // face 0
        indices[0] = 2;     indices[1] = 1;     indices[2] = 3;
        // face 1
        indices[3] = 3;     indices[4] = 1;     indices[5] = 0;
        // LEFT
        // face 2
        indices[6] = 1;     indices[7] = 5;     indices[8] = 0;
        // face 3
        indices[9] = 0;     indices[10] = 5;    indices[11] = 4;
        // BACK
        // face 4
        indices[12] = 1;    indices[13] = 2;    indices[14] = 5;
        // face 5
        indices[15] = 5;    indices[16] = 2;    indices[17] = 6;
        // FRONT
        // face 6
        indices[18] = 3;    indices[19] = 0;    indices[20] = 7;
        // face 7
        indices[21] = 7;    indices[22] = 0;    indices[23] = 4;
        // RIGHT
        // RIGHT 8
        indices[24] = 7;    indices[25] = 2;    indices[26] = 3;
        // face 9
        indices[27] = 2;    indices[28] = 7;    indices[29] = 6;
        // BOTTOM
        // face 10
        indices[30] = 7;    indices[31] = 5;    indices[32] = 6;
        // face 11
        indices[33] = 5;    indices[34] = 7;    indices[35] = 4;
        
        m.build( vertices.size( ), reinterpret_cast< float * >( &vertices[0] ),
                faces.size( ), &faces[0], &indices[0] );
    }
    
    // Creates the quad mesh of a cube
    void create_quads_cube( HMesh::Manifold& m, double edge_length )
    {
        vector< Vec3f > vertices( 8 );
        vector< int >   faces   ( 6, 4 );
        vector< int >   indices ( 24 );
        
        create_cube_vertices(edge_length, vertices);
        
        // for each face save the indices of its vertices
        // UP
        indices[0]  = 2;    indices[1]  = 1;    indices[2]  = 0;    indices[3]  = 3;
        // LEFT
        indices[4]  = 1;    indices[5]  = 5;    indices[6]  = 4;    indices[7]  = 0;
        // BACK
        indices[8]  = 1;    indices[9]  = 2;    indices[10] = 6;    indices[11] = 5;
        // FRONT
        indices[12] = 3;    indices[13] = 0;    indices[14] = 4;    indices[15] = 7;
        // RIGHT
        indices[16] = 7;    indices[17] = 6;    indices[18] = 2;    indices[19] = 3;
        // BOTTOM
        indices[20] = 7;    indices[21] = 4;    indices[22] = 5;    indices[23] = 6;
        
        m.build( vertices.size( ), reinterpret_cast< float * >( &vertices[0] ),
                faces.size( ), &faces[0], &indices[0] );
    }
}
}