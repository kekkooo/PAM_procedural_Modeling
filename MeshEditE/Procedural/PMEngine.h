//
//  PMEngine.h
//  MeshEditE
//
//  Created by Francesco Usai on 30/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__PMEngine__
#define __MeshEditE__PMEngine__

#include <iostream>
#include <GEL/HMesh/Manifold.h>
#include <polarize.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>


using namespace HMesh;
using namespace std;

namespace Procedural
{
    // For now it is able to manage only one single connected component.
    // If I want to have different components I need to navigate the skeleton
    // and extract the different components and keep track of them.
    // each connected component needs to have its polesList and edge_info_container
    
    // This structure keeps the list of all the poles of the mesh for fast indexing
    // after each structural operation on the mesh it must be invalidated and updated
    struct PolesList
    {
    private:
        vector< VertexID >  poles;
        bool                is_valid;
        // operations
    public:
        PolesList() { is_valid = false; }
        
        inline vector< VertexID > const Poles()         { return poles; }
        inline bool               const IsValid()       { return is_valid; }
        inline void                     Invalidate()    { is_valid = false; }
        
        // updates the pole list with the input mesh
        void Update( Manifold *mesh, bool evenIfIsValid = false )
        {
            if( is_valid && !evenIfIsValid ) return;
            poles.clear();
            for( VertexID vid : mesh->vertices() )
            {
                if( is_pole( *mesh, vid )) { poles.push_back( vid ); }
            }
            is_valid = true;
        }
    };
    
    struct EdgeInfoContainer
    {
    private:
        HalfEdgeAttributeVector<EdgeInfo> edge_info;
        bool                              is_valid;
    public:
        EdgeInfoContainer() { is_valid = false; }
        
        inline HalfEdgeAttributeVector<EdgeInfo>    const edgeInfo()    { return edge_info; }
        inline bool                                 const IsValid()     { return is_valid;  }
        inline void                                       Invalidate()  { is_valid = false; }
        
        // updates the pole list with the input mesh
        void Update( Manifold *mesh, bool also_junctions = false, bool evenIfIsValid = false )
        {
            if( is_valid && !evenIfIsValid ) return;
            edge_info = label_PAM_edges( *mesh );
            if( also_junctions )
                Procedural::Structure::LabelJunctions( *mesh, edge_info );
            is_valid = true;
        }
        void UpdateWithJunctions( Manifold *mesh )
        {
            if( !is_valid ) Update( mesh, true );
            else Procedural::Structure::LabelJunctions( *mesh, edge_info );
        }
        
    };
    
    
    
    

class Engine
{
    
    
public :
                        Engine                          ( );
            void        buildCube                       ( );
            void        buildCube                       ( Manifold& mesh);
            void        smooth                          ( Manifold& mesh);
            void        extrudePoles                    ( );
            void        polarSubdivision                ( );
    

            void        setMesh( Manifold* mesh);
    inline  Manifold*   getMesh()                               { return m; }
    inline  void        initMesh( )                             { m = new Manifold(); }

private :
    Manifold            *m;
    PolesList           _polesList;
    EdgeInfoContainer   _edges_info_container;
    
    inline  void        invalidatePolesList() { _polesList.Invalidate(); }
    inline  void        invalidateEdgeInfo()  { _edges_info_container.Invalidate(); }
    inline void         invalidateAll()       { invalidateEdgeInfo(); invalidatePolesList(); }

    
    
};

}

#endif /* defined(__MeshEditE__PMEngine__) */
