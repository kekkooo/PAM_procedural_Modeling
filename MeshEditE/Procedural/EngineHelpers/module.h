//
//  module.h
//  MeshEditE
//
//  Created by Francesco Usai on 24/09/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__module__
#define __MeshEditE__module__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <GEL/CGLA/Vec3d.h>

namespace Procedural{
    namespace EngineHelpers{
        
        typedef HMesh::VertexID             node;
        typedef CGLA::Vec3d                 normal;
        typedef std::pair< node, node >     arc;
        typedef std::map< arc, double >     distance_map;
        typedef std::map< arc, double >     angle_map;

        struct ModuleDescriptor
        {
            // poles and normals are in sync
            std::vector<node>   poles;
            std::vector<normal> normals;
            distance_map        dist;
            angle_map           angles;
            double              max_dist;
        };
        
        struct Module{
            HMesh::Manifold  *m;
            ModuleDescriptor d;
        };
        
        static Module buildModule( HMesh::Manifold& mesh )
        {
            Module m;
            m.m = &mesh;
            m.d.max_dist = numeric_limits<double>::min();
            // detect the poles
            for( node vid : m.m->vertices() )
            {
                if( is_pole( *m.m, vid )) { m.d.poles.push_back( vid ); }
            }
            // calculate the normals
            for( node pole : m.d.poles )
            {
                normal n = Procedural::Geometry::vertex_normal( *m.m, pole );
                n.normalize();
                m.d.normals.push_back(n);
            }
            // calculate distance and angle distance between normals
            for( int i = 0; i < m.d.poles.size(); ++i )
            {
                for( int j = i+1; j < m.d.poles.size(); ++j )
                {
                    arc key     = std::make_pair(m.d.poles[i], m.d.poles[j] );
                    double dist = ( m.m->pos( m.d.poles[i] ) - m.m->pos( m.d.poles[j] )).length();
                    if( dist > m.d.max_dist ) { m.d.max_dist = dist; }
                    m.d.dist  [ key ] = dist;
                    m.d.angles[ key ] = CGLA::dot( m.d.normals[i], m.d.normals[j] );
                }
            }            
            return m;
        }
    }
    
}

#endif /* defined(__MeshEditE__module__) */
