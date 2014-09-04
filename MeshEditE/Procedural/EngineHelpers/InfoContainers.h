//
//  InfoContainers.h
//  MeshEditE
//
//  Created by Francesco Usai on 01/09/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef MeshEditE_InfoContainers_h
#define MeshEditE_InfoContainers_h

#include <GEL/HMesh/Manifold.h>
#include <polarize.h>
#include <MeshEditE/Procedural/Helpers/structural_helpers.h>
#include <MeshEditE/Procedural/Helpers/geometric_properties.h>
#include <unordered_set>


using namespace std;
using namespace HMesh;

namespace Procedural {
    namespace EngineHelpers{
        
typedef VertexAttributeVector< Procedural::Geometry::DistanceMetrics >  DistanceVector;
typedef VertexAttributeVector<double>                                   AngleVector;
typedef HalfEdgeAttributeVector<double>                                 LengthVector;

        
// For now it is able to manage only one single connected component.
// If I want to have different components I need to navigate the skeleton
// and extract the different components and keep track of them.
// each connected component needs to have its polesList and edge_info_container

// This structure keeps the list of all the poles of the mesh for fast indexing
// after each structural operation on the mesh it must be invalidated and updated
struct PolesList
{
private:
    vector< VertexID >                  poles;
    map< VertexID, int >                pole_valency;
    map< VertexID, int >                pole_age;
    bool                                is_valid;
    size_t                              no_poles;
    // operations
    int                         const   indexOfPole( VertexID pole )
    {
        assert( IsPole( pole ));
        auto iter = std::find( poles.begin(), poles.end(), pole );
        return iter - poles.begin();
    }
    
public:
    PolesList() { is_valid = false; no_poles = 0; }
    
    inline vector< VertexID >   const   Poles()         { return poles;     }
    inline bool                 const   IsValid()       { return is_valid;  }
    inline void                         Invalidate()    { is_valid = false; }
    inline size_t               const   No_Poles()      { return no_poles;  }
    inline bool                 const   IsPole( VertexID vertex )
    {
        return std::find( poles.begin(), poles.end(), vertex ) != poles.end();
    }
           int                  const   PoleValency( VertexID pole )
    {
        assert( IsPole( pole ));
        return pole_valency[pole];
    }
    inline int                  const   MeanPoleValency()
    {
        int sum = 0;
        for( auto v : pole_valency ) { sum += v.second; }
        return sum/No_Poles();
    }
    inline int                  const   PoleAge( VertexID pole )
    {
        assert( IsPole(pole));
        return pole_age[pole];
    }
    
    // updates the pole list with the input mesh
    void Update( Manifold *mesh, bool evenIfIsValid = true )
    {
        if( is_valid && !evenIfIsValid ) return;
        for( VertexID vid : mesh->vertices() )
        {
            // if is a pole for the mesh
            if( is_pole( *mesh, vid ))
            {
                // it was a pole in the last iteration
                if( IsPole( vid ))  { pole_age[vid]++;   }
                else
                {
                    poles.push_back( vid );
                    pole_age[vid] = 0;
                }
                pole_valency[vid] = valency( *mesh, vid );
            }
            else
            {
                if( IsPole( vid ))
                {
                    // if it was a pole in the last iteration, but it is not anymore
                    poles.erase( std::find( poles.begin(), poles.end(), vid) );
                    pole_valency.erase( vid );
                    pole_age.erase( vid );
                }
            }
        }
        is_valid = true;
        no_poles = poles.size();
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
    void Update( Manifold *mesh, bool also_junctions = true, bool evenIfIsValid = true )
    {
//        if( is_valid && !evenIfIsValid ) return;
        edge_info = label_PAM_edges( *mesh );
//        if( also_junctions )
            Procedural::Structure::LabelJunctions( *mesh, edge_info );
        is_valid = true;
    }
    void UpdateWithJunctions( Manifold *mesh )
    {
        if( !is_valid ) Update( mesh, true );
        Procedural::Structure::LabelJunctions( *mesh, edge_info );
    }        
};
        
struct GeometricInfoContainer
{
private:
    DistanceVector      pole_dist;
    DistanceVector      junction_dist;
    DistanceVector      combined_dist;
    AngleVector         angles;
    AngleVector         spine_angles;
    AngleVector         rib_angles;
    LengthVector        edge_lengths;
    double              mean_length;    
    bool                is_valid;
public:
    inline DistanceVector   const PoleDistance()       { return pole_dist;     }
    inline DistanceVector   const JunctionDistance()   { return junction_dist; }
    inline DistanceVector   const CombinedDistance()   { return combined_dist; }
    inline AngleVector      const CombinedAngles()     { return angles;        }
    inline AngleVector      const RibAngles()          { return rib_angles;    }
    inline AngleVector      const SpineAngles()        { return spine_angles;  }
    inline LengthVector     const EdgeLength()         { return edge_lengths;  }
    inline double           const MeanLength()         { return mean_length;   }
    
    inline bool             const IsValid()             { return is_valid;  }
    inline void                   Invalidate()          { is_valid = false; }
    
    // maybe the right thing to do is to have a function for each group, in order to
    // update them separately, because they are not always needed and are expensive
    // also I need a IsValid for each container.
    // probably a good Idea to create a container for each geometric info
    void Update( Manifold *mesh, EdgeInfoContainer &edgeInfo,  bool evenIfIsValid = true )
    {
        if( is_valid && !evenIfIsValid ) return;

        Procedural::Geometry::dihedral_angles                   ( *mesh, edgeInfo.edgeInfo(), angles                 );
        Procedural::Geometry::dihedral_angles                   ( *mesh, edgeInfo.edgeInfo(), spine_angles,  SPINE   );
        Procedural::Geometry::dihedral_angles                   ( *mesh, edgeInfo.edgeInfo(), rib_angles,    RIB     );

        Procedural::Geometry::distance_from_poles               ( *mesh, edgeInfo.edgeInfo(), pole_dist,     false   );
        Procedural::Geometry::distance_from_junctions           ( *mesh, edgeInfo.edgeInfo(), junction_dist, false   );
        Procedural::Geometry::distance_from_poles_and_junctions ( *mesh, edgeInfo.edgeInfo(), combined_dist, false   );
        UpdateEdgeLengths ( *mesh );
        is_valid = true;
    }
    
    void UpdateEdgeLengths( Manifold& mesh )
    {
        mean_length =
        Procedural::Geometry::edge_length                       ( mesh, edge_lengths                               );
    }
};
        
        
struct PoleTracking
{
private:
        map< VertexID, VertexID > matches;
    
public:
    inline bool         ExistsMatching( VertexID v ) { return matches.count(v) > 0; }
    inline void         Clear( )                     { matches.clear(); }
           VertexID     GetMatch( VertexID v )
    {
        if( !ExistsMatching(v)) return InvalidVertexID;
        return matches[v];
    }
    inline void        AddMatch( VertexID v, VertexID match ) { matches[v] = match; }
};
        
        
struct PoleTrajectory
{
    int         no_calls;
    int d1;
    int d2;
    CGLA::Vec3d current_dir;
    double      total_length;
    PoleTrajectory() { no_calls = 0; total_length = 0.0; }
};
        
        
}}


#endif
