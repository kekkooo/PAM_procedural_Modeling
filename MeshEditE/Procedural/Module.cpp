//
//  Module.cpp
//  MeshEditE
//
//  Created by Francesco Usai on 24/04/15.
//  Copyright (c) 2015 J. Andreas Bærentzen. All rights reserved.
//

#include "Module.h"
#include "Test.h"
#include <fstream>
#include <streambuf>
#include <iostream>

#include "rapidjson/document.h"

#include <GEL/HMesh/obj_load.h>

#include "MeshEditE/Procedural/Helpers/geometric_properties.h"

#include "Plane.h"
#include "eigenv.h"
#include "Test.h"

using namespace std;
using namespace HMesh;
using namespace CGLA;
using namespace Procedural::Geometry;

namespace Procedural{

Module::Module( std::string path, std::string config, Moduletype mType ){
    cout << "trying to load file : " << path << endl;
    ifstream f( path );
    assert( f.good() );
    
    this->m = new Manifold();
    m->clear();
    
    obj_load( path, *this->m );
    bsphere( *m, bsphere_center, bsphere_radius );
    this->type = mType;
    BuildPoleInfo();
    LoadPoleConfig( config );
    
    this->skeleton = new Skeleton();
    this->skeleton->build( *m, this->poleSet );
    this->skeleton->saveToFile("//Users//francescousai//Desktop//example.skel");
}
    
Module::Module( Manifold &manifold, Moduletype mType ){
    this->m = &manifold;
    Vec3d centroid;
    double radius;
    bsphere( *this->m, centroid, radius );
    this->bsphere_center = centroid;
    this->bsphere_radius = radius;
    
    BuildPoleInfo();
    
    this->skeleton = new Skeleton();
    this->skeleton->build( *m, this->poleSet );
}

void Module::LoadPoleConfig( std::string path ){
    std::cout << path;
    std::ifstream t( path );
//    assert( t.good() );
    if( !t.std::__1::ios_base::good()){
        cout << "configuration file " << path << "not found or not good" << endl;
        return;
    }
    std::string json( (std::istreambuf_iterator<char>(t)),
                     std::istreambuf_iterator<char>());
    
    rapidjson::Document d;
    d.Parse( json.c_str() );
    rapidjson::Value &tb = d["config"];
    assert( tb.IsArray() );
    
    for( rapidjson::SizeType i = 0; i < tb.Size(); ++i ){
        assert( tb[i].HasMember( "pole_id" ));
        assert( tb[i].HasMember( "direction" )); // pole neighbor that expresses the preferred direction of alignment
        assert( tb[i].HasMember( "bilateral" ));
        assert( tb[i].HasMember( "active" ));
        assert( tb[i].HasMember( "connect_to_self" ));
        
        int     pole_id         = tb[i]["pole_id"].GetInt();
        int     direction       = tb[i]["direction"].GetInt();
        bool    bilateral       = tb[i]["bilateral"].GetBool();
        bool    active          = tb[i]["active"].GetBool();
        bool    connect_to_self = tb[i]["connect_to_self"].GetBool();
        
        bool done = false;
        for( VertexID pole : poleList ){
            if( done ){ assert(pole_id != pole.get_index()); continue; }
            if( pole_id == pole.get_index()){
                poleInfoMap[pole].can_connect_to_self           = connect_to_self;
                poleInfoMap[pole].isActive                      = active;
                poleInfoMap[pole].anisotropy.is_bilateral       = bilateral;
                
                if( direction >= 0 ){
                    poleInfoMap[pole].anisotropy.is_defined    = true;
                    
                    VertexID neighbor = InvalidVertexID;
                    // find which neighbor is the one
                    for( Walker w = m->walker( pole ); (( !w.full_circle( )) && ( neighbor == InvalidVertexID )) ; w = w.circulate_vertex_ccw() ){
                        if( w.vertex().get_index() == direction ){ neighbor = w.vertex(); }
                    }
                    assert( neighbor != InvalidVertexID );

                    poleInfoMap[pole].anisotropy.directionID    = neighbor;
                    getPoleAnisotropy( poleInfoMap[pole].geometry.pos, poleInfoMap[pole].geometry.normal,
                                      m->pos( poleInfoMap[pole].anisotropy.directionID ),
                                      poleInfoMap[pole].anisotropy.direction );


                }
                done = true;
            }
        }
        assert( done ); // means that a vertex with id equal to pole_id was found
        
    }

}
    
void Module::BuildPoleInfo(){
    
    assert( m != NULL );
    
    for( VertexID vid : m->vertices() ){
        if( is_pole( *m, vid )){
            PoleInfo pi;
            pi.moduleType       = this->type;
            pi.original_id      = vid;
            pi.geometry.valence = valence( *m, vid );
            pi.geometry.pos     = m->pos( vid );
            Vec3d n             = vertex_normal( *m, vid );
            n.normalize();
            truncateVec3d( n );
            pi.geometry.normal  = n;

            this->poleList.push_back( vid );
            this->poleSet.insert( vid );
            this->poleInfoMap[vid] = pi;
        }
    }
}
    

void Module::getPoleAnisotropy( const CGLA::Vec3d& pole, const CGLA::Vec3d& normal, const CGLA::Vec3d& neighbor, CGLA::Vec3d& dir ){    
    
    Plane p( pole, normal );
    dir = p.projectDirection( pole, neighbor );
    truncateVec3d( dir );
}

    
void Module::updateDirections( const HMesh::Manifold& main_mesh ){
    
    for( VertexID v : this->poleList ){
        CGLA::Vec3d normal = vertex_normal( main_mesh, v );
        normal.normalize();
        poleInfoMap[v].geometry.normal = normal ;
        truncateVec3d( poleInfoMap[v].geometry.normal );
        poleInfoMap[v].geometry.pos = main_mesh.pos( v ) ;
        
        getPoleAnisotropy( poleInfoMap[v].geometry.pos, poleInfoMap[v].geometry.normal,
                           main_mesh.pos( poleInfoMap[v].anisotropy.directionID ),
                           poleInfoMap[v].anisotropy.direction );

    }        
}
    
Module& Module::getTransformedModule( const CGLA::Mat4x4d &T, bool transform_geometry )
{
    Module *M = new Module();
    M->m             = this->m;
    M->poleList      = PoleList();
    M->poleInfoMap   = PoleInfoMap();
    
    M->no_of_glueings = this->no_of_glueings;
    M->bsphere_center = T.mul_3D_point( bsphere_center );
    M->bsphere_radius = bsphere_radius;
    
    for( VertexID vid : this->poleList ){
        M->poleList.push_back( vid );
        M->poleInfoMap[vid].original_id             = poleInfoMap[vid].original_id;
        M->poleInfoMap[vid].geometry.pos            = T.mul_3D_point( poleInfoMap[vid].geometry.pos );
        M->poleInfoMap[vid].geometry.normal         = mul_3D_dir( T, poleInfoMap[vid].geometry.normal );
        M->poleInfoMap[vid].geometry.normal.normalize();
        M->poleInfoMap[vid].geometry.valence        = poleInfoMap[vid].geometry.valence;

        M->poleInfoMap[vid].anisotropy.is_defined   = poleInfoMap[vid].anisotropy.is_defined;
        M->poleInfoMap[vid].anisotropy.direction    = mul_3D_dir( T, poleInfoMap[vid].anisotropy.direction );
        M->poleInfoMap[vid].anisotropy.direction.normalize();
        
        M->poleInfoMap[vid].anisotropy.is_bilateral = poleInfoMap[vid].anisotropy.is_bilateral;
        M->poleInfoMap[vid].anisotropy.directionID  = poleInfoMap[vid].anisotropy.directionID;
        
        M->poleInfoMap[vid].moduleType              = poleInfoMap[vid].moduleType;
        M->poleInfoMap[vid].age                     = poleInfoMap[vid].age;
        M->poleInfoMap[vid].isFree                  = poleInfoMap[vid].isFree;
        M->poleInfoMap[vid].can_connect_to_self     = poleInfoMap[vid].can_connect_to_self;
        M->poleInfoMap[vid].isActive                = poleInfoMap[vid].isActive;
        
        // truncate directions
        truncateVec3d( M->poleInfoMap[vid].geometry.normal );
        truncateVec3d( M->poleInfoMap[vid].anisotropy.direction );
        
        bool assert_pos =
            isnan( M->poleInfoMap[vid].geometry.pos[0] ) ||
            isnan( M->poleInfoMap[vid].geometry.pos[1] ) ||
            isnan( M->poleInfoMap[vid].geometry.pos[2] );
        
        bool assert_normal =
            isnan( M->poleInfoMap[vid].geometry.normal[0] ) ||
            isnan( M->poleInfoMap[vid].geometry.normal[1] ) ||
            isnan( M->poleInfoMap[vid].geometry.normal[2] );
        
        if( M->poleInfoMap[vid].anisotropy.is_defined ){
            bool assert_anis =
                isnan( M->poleInfoMap[vid].anisotropy.direction[0] ) ||
                isnan( M->poleInfoMap[vid].anisotropy.direction[1] ) ||
                isnan( M->poleInfoMap[vid].anisotropy.direction[2] );
            assert( !assert_anis );
        }
        
        assert( !assert_pos );
        assert( !assert_normal );
    }
    
    if( transform_geometry ){
        for( VertexID v : m->vertices()){
            m->pos( v ) = T.mul_3D_point( m->pos( v ));
        }
    }
    
    M->sanityCheck();
    
    assert( M->poleInfoMap.size() == this->poleInfoMap.size( ));
    
    M->skeleton = new Skeleton();    
    M->skeleton->copyNew( *skeleton );
    M->skeleton->transform( T );
    
    return *M;
}
    
void Module::reAlignIDs(HMesh::VertexIDRemap &remapper){
    // realign poleList
    PoleInfoMap p;
    for( int i = 0; i < poleList.size(); ++i ){
        VertexID oldID = poleList[i];
        VertexID aniID = poleInfoMap[oldID].anisotropy.directionID;
        
        VertexID newAnisID = remapper[aniID];
        VertexID newID      = remapper[oldID];
        
        assert( newAnisID != InvalidVertexID );
        assert( newID != InvalidVertexID );
        
        p[newID]       = poleInfoMap[poleList[i]];
        poleList[i]    = newID;
        p[newID].anisotropy.directionID = newAnisID;
    }
    poleInfoMap = std::move( p );
    Skeleton* temp = skeleton;
    skeleton = new Skeleton();
    skeleton->copyAndRealignIDs( *temp, remapper );
}

const PoleInfo& Module::getPoleInfo( HMesh::VertexID p ) const{
    assert(poleInfoMap.count( p ) > 0);
    return poleInfoMap.at(p);
}
    
bool Module::isPole( HMesh::VertexID v ){
    return (( find( poleList.begin(), poleList.end(), v ) != poleList.end()) && ( poleInfoMap.count(v) > 0 ) );
}
    
const Skeleton& Module::getSkeleton() const{
    return *( this->skeleton );
}
    
/*** STATIC ***/
// this should be elsewhere, but I don't know where to put it.
bool Module::poleCanMatch( const PoleInfo& p1, const PoleInfo& p2){
//    cout << "testing : " << p1.moduleType << " -> " << p1.original_id << " with "
//                         << p2.moduleType << " -> " << p2.original_id << endl;

    if( !p1.isActive ){
//        cout << "KO ( " << p1.moduleType << ", " << p1.original_id << ") is not active" << endl;
        return false;
    }
    
    if( !p2.isActive ){
//        cout << "KO ( " << p2.moduleType << ", " << p2.original_id << ") is not active" << endl;
        return false;
    }

    
    // test valence
    if( p1.geometry.valence != p2.geometry.valence ) {
//        cout << "KO for valence (" << p1.geometry.valence << ", " << p2.geometry.valence << " )" << endl;
        return false;
    }

    // test self connectability
    if( p1.moduleType == p2.moduleType ){
        if( p1.original_id == p2.original_id ){
            if( !( p1.can_connect_to_self && p2.can_connect_to_self )){
//                cout << "KO connect_to_self" << endl;
                return false;
            }
        }
    }
    return true;
}
    
    void Module::sanityCheck(){
        for( const auto& item : poleInfoMap ){
            checkVec3( item.second.geometry.normal );
            if( item.second.anisotropy.is_defined){ checkVec3( item.second.anisotropy.direction ); }
        }
    }
    

}