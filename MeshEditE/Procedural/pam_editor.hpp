//
//  pam_editor.hpp
//  MeshEditE
//
//  Created by Francesco Usai on 10/02/16.
//  Copyright © 2016 J. Andreas Bærentzen. All rights reserved.
//

#ifndef pam_editor_hpp
#define pam_editor_hpp

#include <stdio.h>
#include <set>
#include <map>
#include <GEL/GLGraphics/ManifoldRenderer.h>
#include <GEL/HMesh/Manifold.h>
#include <GEL/CGLA/Vec3d.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>

#include <polarize.h>
#include <Module.h>
#include <test.h>


namespace Procedural{
    
    typedef CGAL::Simple_cartesian<double>                                                  K;
    typedef K::Point_3                                                                      CGAL_Point;
    typedef K::Triangle_3                                                                   CGAL_Triangle;
    typedef K::Line_3                                                                       CGAL_Line;
    typedef K::Segment_3                                                                    CGAL_Segment;
    typedef K::Plane_3                                                                      CGAL_Plane;
    typedef K::Vector_3                                                                     CGAL_Vector;
    typedef K::Direction_3                                                                  CGAL_Direction;

    
enum cross_section_type{ ROUND, UNILATERAL, BILATERAL };
    
    struct pole_plane_info{
        CGAL_Point      pole_pos;
        CGAL_Vector     pole_direction;
        CGAL_Plane      pole_plane;
        std::map< HMesh::VertexID, CGLA::Vec3d > pole_neighbourhood;
        std::map< HMesh::VertexID, CGLA::Vec3d > pole_neighbourhood_vectors;
        CGLA::Mat4x4d rot_cw;
        CGLA::Mat4x4d rot_ccw;

    };
    
class PamEditor{
public:
    // singleton implementation
    static  PamEditor& getCurrentPamEditor(){
        static PamEditor instance;
        return instance;
    }

    PamEditor( PamEditor const& )           = delete;
    void operator   = (PamEditor const&)    = delete;
    PamEditor(){ }
    //
    
    
    Procedural::Module* module = nullptr;
    HMesh::VertexID             single_selected_pole = HMesh::InvalidVertexID;
    std::set<HMesh::VertexID>   selected_poles;
    std::map<HMesh::VertexID, pole_plane_info> pole_to_plane_info;
    

    inline void setCrossSection( cross_section_type cst ){
        if( single_selected_pole == HMesh::InvalidVertexID ){
            std::cout << "not rotating. select exactly 1 pole" << std::endl;
            return;
        }
        pole_plane_info& info = pole_to_plane_info[single_selected_pole];
        
        switch( cst ){
            case  ROUND : {
                
                module->poleInfoMap[single_selected_pole].anisotropy.is_defined     = false;
                module->poleInfoMap[single_selected_pole].anisotropy.is_bilateral   = false;
                break;
            }
            case UNILATERAL : {
                module->poleInfoMap[single_selected_pole].anisotropy.is_defined     = true;
                module->poleInfoMap[single_selected_pole].anisotropy.is_bilateral   = false;
                module->poleInfoMap[single_selected_pole].anisotropy.directionID    = (*info.pole_neighbourhood_vectors.begin()).first;
                module->poleInfoMap[single_selected_pole].anisotropy.direction      = (*info.pole_neighbourhood_vectors.begin()).second;
                break;
            }
            case BILATERAL : {
                module->poleInfoMap[single_selected_pole].anisotropy.is_defined     = true;
                module->poleInfoMap[single_selected_pole].anisotropy.is_bilateral   = true;
                module->poleInfoMap[single_selected_pole].anisotropy.directionID    = (*info.pole_neighbourhood_vectors.begin()).first;
                module->poleInfoMap[single_selected_pole].anisotropy.direction      = (*info.pole_neighbourhood_vectors.begin()).second;
                break;
            }
            default:
                std::cout << "there should not be other options to this " << std::endl;
                assert(false);
        }
        
    }
    
    // this is going to be called after the module.pole.anisotropy.direction has been changed
    inline void updateVerticesAfterRotation( CGLA::Mat4x4d& T){
        assert( single_selected_pole != HMesh::InvalidVertexID );
        double              angle_dist      = std::numeric_limits<double>::max();
        HMesh::VertexID     closest         = HMesh::InvalidVertexID;
        pole_plane_info&    info            = pole_to_plane_info[single_selected_pole];
        
        CGLA::Vec3d         curr_neighbor   = module->m->pos( module->getPoleInfo( single_selected_pole ).anisotropy.directionID );
        CGAL_Point          cgal_neighbor( curr_neighbor[0], curr_neighbor[1], curr_neighbor[2] );
        CGAL_Point          projected       = info.pole_plane.projection( cgal_neighbor );
        CGAL_Vector         cgal_vector( info.pole_pos, projected );
        CGLA::Vec3d         curr_direction( cgal_vector[0], cgal_vector[2], cgal_vector[2] );

        
        for( auto& item : info.pole_neighbourhood_vectors ){
            double cosine = dot( curr_direction, item.second );
            if(( 1.0 - cosine ) < angle_dist ){
                angle_dist = 1.0 - cosine;
                closest = item.first;
            }
        }
        assert( closest != HMesh::InvalidVertexID );
        if( closest == module->getPoleInfo( single_selected_pole ).anisotropy.directionID ){
            // in this case I should rotate the points to accomodate the new direction
            // update the stored directions
            for( auto& item : info.pole_neighbourhood ){
                module->m->pos( item.first ) = T.mul_3D_point( module->m->pos( item.first ));
            }
            
        }else{
            // in this case I should change DirectionID
            module->poleInfoMap[single_selected_pole].anisotropy.directionID = closest;
            // rollback the neighborhood to the original positions
            for( auto& item : info.pole_neighbourhood){
                module->m->pos( item.first ) = item.second;
            }
            // compute the new configuration
        }
        
        // in each case, I must update directions
        for( auto& item : info.pole_neighbourhood ){
            CGLA::Vec3d neighbor = module->m->pos( item.first );
            CGAL_Point  cgal_neighbor( neighbor[0], neighbor[1], neighbor[2] );
            CGAL_Point  projected = info.pole_plane.projection( cgal_neighbor );
            CGAL_Vector plane_vector( info.pole_pos, projected );
            info.pole_neighbourhood_vectors[item.first] = CGLA::Vec3d( plane_vector[0], plane_vector[1], plane_vector[2] );
        }        
    }
    
    // le rotazioni sono sbagliate perché devo ruotare la direzione,
    // e in seguito devo ruotare i vertici per adeguarli alla nuova posizione
    
    inline void rotate_selected_cw(){
        if( single_selected_pole == HMesh::InvalidVertexID ){
            std::cout << "not rotating. select only 1 pole" << std::endl;
            return;
        }
        if( !module->poleInfoMap.at( single_selected_pole ).anisotropy.is_defined ){ return; }

        pole_plane_info& info = pole_to_plane_info[single_selected_pole];
        module->poleInfoMap[single_selected_pole].anisotropy.direction =
            info.rot_cw.mul_3D_vector( module->poleInfoMap[single_selected_pole].anisotropy.direction );
        
        updateVerticesAfterRotation( info.rot_cw );
        debugColorization();
    }
        
    inline void rotate_selected_ccw(){
        if( single_selected_pole == HMesh::InvalidVertexID ){
            std::cout << "not rotating. select only 1 pole" << std::endl;
            return;
        }
        pole_plane_info& info = pole_to_plane_info[single_selected_pole];
        module->poleInfoMap[single_selected_pole].anisotropy.direction =
            info.rot_ccw.mul_3D_vector( module->poleInfoMap[single_selected_pole].anisotropy.direction );
        
        updateVerticesAfterRotation( info.rot_ccw );
        debugColorization();
    }

    // QUA DEVO PROIETTARE I PUNTI SUL PIANO
    inline void computePlanesInfo(){
        for( HMesh::VertexID vid : module->poleSet ){
            CGLA::Vec3d pos     = module->getPoleInfo( vid ).geometry.pos;
            CGLA::Vec3d normal  = module->getPoleInfo( vid ).geometry.normal;
            
            CGAL_Point  p( pos[0],    pos[1],    pos[2]    );
            CGAL_Vector n( normal[0], normal[1], normal[2] );
            CGAL_Plane  plane( p, n );
            pole_to_plane_info[vid].pole_pos        = p;
            pole_to_plane_info[vid].pole_direction  = n;
            pole_to_plane_info[vid].pole_plane      = plane;
            
            CGLA::Mat4x4d to_origin     = CGLA::translation_Mat4x4d( -pos );
            CGLA::Mat4x4d back_to_pos   = CGLA::translation_Mat4x4d( pos );
            
            pole_to_plane_info[vid].rot_cw          = back_to_pos * get_rotation_mat4d( normal,  0.174533 ) * to_origin;
            pole_to_plane_info[vid].rot_ccw         = back_to_pos * get_rotation_mat4d( normal, -0.174533 ) * to_origin;

            HMesh::Walker w = module->m->walker( vid );


            for( ; !w.full_circle(); w = w.circulate_vertex_ccw() ){
                CGLA::Vec3d neighbor     = module->m->pos( w.vertex( ));
                CGLA::Vec3d neighbor_dir = neighbor - pos;
                pole_to_plane_info[vid].pole_neighbourhood[w.vertex()]         = neighbor;
                pole_to_plane_info[vid].pole_neighbourhood_vectors[w.vertex()] = neighbor_dir;
            }
        }
    }
    
    inline void set_mesh( HMesh::Manifold* mesh ){
        if( module != nullptr ){ delete( module ); }
        
        module = new Module();
        module->m = mesh;
        module->BuildPoleInfo();
    }
    
    inline void set_config( std::string path){
        module->LoadPoleConfig( path );
        computePlanesInfo();
        debugColorization();
    }
    
    inline void clear_selection(){
        selected_poles.clear();
        single_selected_pole = HMesh::InvalidVertexID;
    }
    
    inline void add_to_selection( HMesh::VertexID pole ){
        assert( module    != nullptr );
        assert( module->m != nullptr );
        
        if( module->poleSet.count(pole) > 0 ){
            selected_poles.insert( pole );
            if( selected_poles.size() ==  1 ){
                single_selected_pole = pole;
            }else{
                single_selected_pole = HMesh::InvalidVertexID;
            }
            
        }
    }

    inline void debugColorization(){
        
        CGLA::Vec3f v_color( 0.0, 0.0, 1.0 );
        CGLA::Vec3f f_color( 0.9, 0.9, 0.9 );
        CGLA::Vec3f e_color( 0.0, 0.0, 0.0 );
        CGLA::Vec3f d_color( 1.0, 0.0, 0.0 ); // directions
        CGLA::Vec3f p_color( 0.0, 1.0, 1.0 ); // pole
        for( auto vit = module->m->vertices().begin(); vit != module->m->vertices().end(); ++vit )
        {
            assert( module->m->in_use( *vit ));
            bool is_pole = module->isPole( *vit ) > 0;
            
            
            GLGraphics::DebugRenderer::vertex_colors[*vit] = is_pole ? p_color : v_color;
            
            HMesh::Walker w = module->m->walker(*vit);
            for(; !w.full_circle(); w = w.circulate_vertex_ccw())
            {
                if( is_pole && module->getPoleInfo(*vit).anisotropy.is_defined ){
                    if( w.vertex() == module->getPoleInfo(*vit).anisotropy.directionID ){
                        GLGraphics::DebugRenderer::edge_colors[w.halfedge()]        = d_color;
                        GLGraphics::DebugRenderer::edge_colors[w.opp().halfedge()]  = d_color;
                    }
                    else{
                        GLGraphics::DebugRenderer::edge_colors[w.halfedge()]        = e_color;
                        GLGraphics::DebugRenderer::edge_colors[w.opp().halfedge()]  = e_color;
                    }
                }
            }
        }
        for( auto fit = module->m->faces().begin(); fit != module->m->faces().end(); ++fit ){
            GLGraphics::DebugRenderer::face_colors[*fit] = f_color;
        }
    }
    
};
}

#endif /* pam_editor_hpp */
