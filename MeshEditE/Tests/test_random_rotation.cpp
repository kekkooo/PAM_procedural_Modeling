//
//  test_random_rotation.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 11/03/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "test_random_rotation.h"
#include "MeshEditE/Procedural/Helpers/module_alignment.h"
#include "MeshEditE/Procedural/Helpers/manifold_copy.h"
#include "MeshEditE/Procedural/Helpers/geometric_properties.h"

using namespace Tests;
using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Procedural::Helpers;
using namespace Procedural::Geometry;

RandomRotationTest& RandomRotationTest::getInstance()
{
    static RandomRotationTest    instance; // Guaranteed to be destroyed.
    // Instantiated on first use.
    return instance;
}

RandomRotationTest::RandomRotationTest()
{
    this->m = NULL;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    randomizer.seed( seed );
}

void RandomRotationTest::setHost( HMesh::Manifold &manifold )
{
    // set the mesh
    this->m = &manifold;
    // save its vertex IDs
//    for( VertexID vid : this->m->vertices( )) {  this->H_vertices.push_back( vid ); }
    
}
void RandomRotationTest::setModule( HMesh::Manifold &manifold )
{
    VertexSet host_ids, module_ids;
    add_manifold( (*this->m), manifold, H_vertices, M_vertices );
//    H_vertices.clear();
//    M_vertices.clear();
//    std::copy( host_ids.begin(),   host_ids.end(),   std::back_inserter( H_vertices ));
//    std::copy( module_ids.begin(), module_ids.end(), std::back_inserter( M_vertices ));
    for( VertexID vid : this->M_vertices )
    {
        this->M_original_pos_map[ vid ] = this->m->pos( vid );
    }
}

void RandomRotationTest::releaseMesh()
{
    this->m = NULL;
}

void RandomRotationTest::rotate( )
{
    // must be initialized
    float rand_max =  static_cast<float>( randomizer.max( )); //RAND_MAX / ( M_PI * 2.0 );
    float x1 = static_cast<float>( randomizer( )) / rand_max,
          x2 = static_cast<float>( randomizer( )) / rand_max,
          x3 = static_cast<float>( randomizer( )) / rand_max;
    
    CGLA::Mat4x4d rot = ModuleAlignment::random_rotation_matrix_arvo( x1, x2, x3 );
    CGLA::Vec3d tr_dir( Vec3d( x1, x2, x3 ));
    tr_dir.normalize();
    
    CGLA::Mat4x4d tr = CGLA::translation_Mat4x4d( tr_dir );
    CGLA::Mat4x4d T = tr * rot;
    
    for( VertexID vid : this->M_vertices )
    {
        this->m->pos(vid) = T.mul_3D_point( this->m->pos( vid ));
    }
}

void RandomRotationTest::translateOutsideBSphere()
{
    // calculate bounding sphere of the two manifolds
    Vec3d  host_centroid, module_centroid;
    double host_radius, module_radius;
    bsphere( *m, H_vertices, module_centroid,   module_radius );
    bsphere( *m, M_vertices, host_centroid,     host_radius   );
    Vec3d   dir     = module_centroid - host_centroid;
    double  length  = host_radius + module_radius - dir.length();
    dir.normalize();
    dir *= length;
    for( VertexID vid : this->M_vertices )
    {
        this->m->pos(vid) += dir;
    }
}

void RandomRotationTest::resetOriginalModulePosition()
{
    for( VertexID vid : this->M_vertices )
    {
        assert( this->M_original_pos_map.count(vid) > 0 );
        this->m->pos( vid ) = this->M_original_pos_map[ vid ];
    }
}
