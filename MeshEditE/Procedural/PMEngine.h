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
#include <MeshEditE/Procedural/EngineHelpers/InfoContainers.h>
#include <MeshEditE/Procedural/Operations/geometric_operations.h>


using namespace HMesh;
using namespace std;
using namespace Procedural::EngineHelpers;

namespace Procedural {

class Engine
{

public :
                        Engine                          ( );
            void        buildCube                       ( );
            void        buildCube                       ( Manifold& mesh);
            void        smooth                          ( Manifold& mesh);
            void        extrudePoles                    ( );
            void        polarSubdivision                ( );
            void        perturbate                      ( );
            void        smooth_near_junctions           ();
            void        add_module                      ();

    
// THE FUNCTIONS ENCLOSED INTO BETWEEN THOSE TWO COMMENTS MUST REPLACED
// THESE THINGS WILL BE PART OF THE FRAMEWORK BUT THEY SHOULD BE PROCEDURAL NOT RANDOM
void addRandomBranches();
void pickABranchAndScaleIt( int mode );
inline int  no_scaling_modes() { return 4; }

// JUST FOR TESTING PURPOSES
    

            void        setMesh( Manifold* mesh);
    inline  Manifold*   getMesh()                               { return m; }
    inline  void        initMesh( )                             { m = new Manifold(); }

private :
    Manifold                *m;
    PolesList               _polesList;
    EdgeInfoContainer       _edges_info_container;
    GeometricInfoContainer  _geometric_info;
    PoleTracking            _pole_tracking;
    int                     _timestamp;
    vertices_info           _v_info;

    HMesh::VertexAttributeVector<int> vertex_selection;
    
    // maybe it will not be needed anymore
    map< VertexID, PoleTrajectory > trajectories;
    //

    
    inline  void        invalidatePolesList()       { _polesList.Invalidate(); }
    inline  void        invalidateEdgeInfo()        { _edges_info_container.Invalidate(); }
    inline  void        invalidateGeometricInfo()   { _geometric_info.Invalidate(); }
            void        invalidateAll();
            void        buildCleanSelection();
    

    
    
};

}

#endif /* defined(__MeshEditE__PMEngine__) */
