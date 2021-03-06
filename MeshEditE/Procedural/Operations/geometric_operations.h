//
//  geometric_operations.h
//  MeshEditE
//
//  Created by Francesco Usai on 08/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__geometric_operations__
#define __MeshEditE__geometric_operations__

#include <iostream>
#include <vector>
#include <tuple>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>
#include "polarize.h"

using namespace std;


namespace Procedural{
    namespace Operations{
        namespace Geometric{


typedef std::vector< std::tuple< HMesh::VertexID, double, CGLA::Vec3d > >  DistanceDirectionVector;

enum VertexType
{
    REGULAR     = 0x01,
    SINGULAR    = 0x02,
    POLE        = 0x04
};

/**********************************************************************
 *******                    SPLITTING                           *******
**********************************************************************/

vector< HMesh::VertexID >   split_ring_of_quads     ( HMesh::Manifold& m, HMesh::HalfEdgeID h   );
void                        split_ring_of_quads     ( HMesh::Manifold& m, HMesh::HalfEdgeID h, int slices );
// h is a rib edge
vector< HMesh::VertexID >   split_from_pole_to_pole ( HMesh::Manifold& m, HMesh::HalfEdgeID h   );

/**********************************************************************
 *******            RIB EDGE LOOP RADIUS OPERATIONS             *******
 **********************************************************************/


void scale_ring_radius  ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                      double ratio, bool smoothing = false );
void set_ring_radius    ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                      double radius );
void set_ring_radius    ( HMesh::Manifold& m, HMesh::HalfEdgeID h,
                                                      std::vector<double> radii );
void set_ring_radius    ( HMesh::Manifold& m, CGLA::Vec3d barycenter,
                                                      DistanceDirectionVector ddv );
            
/**********************************************************************
 *******                POLE GEOMETRIC OPERATIONS               *******
 **********************************************************************/


void flatten_pole ( HMesh::Manifold& m, HMesh::VertexID pole );
void extrude_pole ( HMesh::Manifold& m, HMesh::VertexID v, double length = -1.0,
                   bool pole_add_quads = false, double radius_preserving_ratio = 1.0 );
void smooth_pole  ( HMesh::Manifold& m, HMesh::VertexID pole,
                    HMesh::HalfEdgeAttributeVector<EdgeInfo> edge_info );

void extrude_pole ( HMesh::Manifold&m, HMesh::VertexID v, CGLA::Vec3d direction,
                    bool add_quads, double scaling );
            
void add_ring_around_pole    ( HMesh::Manifold& m, HMesh::VertexID pole, double scaling );
            
/**********************************************************************
 *******                VERTEX GEOMETRIC OPERATIONS               *******
 **********************************************************************/


void move_vertex             ( HMesh::Manifold& m, HMesh::VertexID v, CGLA::Vec3d dir );

void add_noise               ( HMesh::Manifold& m, int vertex_type_flags,
                                double ratiom, double cutoff );
void add_noise               ( HMesh::Manifold& m, int vertex_type_flags, double ratiom,
                               double cutoff, vector< HMesh::VertexID > vertices );

void add_perpendicular_noise ( HMesh::Manifold& m, vector< HMesh::VertexID > &vertices,
                               double amplitude, double avg_edge_length );

void add_perpendicular_noise ( HMesh::Manifold& m, vector< HMesh::VertexID > &vertices,
                               vector< double > &per_vertex_amplitude, double avg_edge_length );


            
}}}
#endif /* defined(__MeshEditE__geometric_operations__) */
