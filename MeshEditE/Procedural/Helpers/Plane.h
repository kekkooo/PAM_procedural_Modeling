//
//  Plane.h
//  MeshEditE
//
//  Created by Francesco Usai on 04/08/14.
//  Copyright (c) 2014 J. Andreas Bærentzen. All rights reserved.
//

#ifndef __MeshEditE__Plane__
#define __MeshEditE__Plane__

#include <iostream>
#include <GEL/CGLA/Vec3d.h>

#define PLANE_EPSILON 0.00001

using namespace CGLA;

namespace Procedural{
    namespace Geometry{

class Plane
{
    
public:
    Plane(Vec3d point, Vec3d normal)
    {
        normal.normalize();
        
        _a = normal[0];
        _b = normal[1];
        _c = normal[2];
        _d = normal[0] * point[0] +
        normal[1] * point[1] +
        normal[2] * point[2];
    }
    
    
    Plane(Vec3d p0, Vec3d p1, Vec3d p2)
    {
        double dx1 = p1[0] - p0[0];
        double dx2 = p2[0] - p0[0];
        double dy1 = p1[1] - p0[1];
        double dy2 = p2[1] - p0[1];
        double dz1 = p1[2] - p0[2];
        double dz2 = p2[2] - p0[2];
        
        double A = (dy1 * dz2) - (dy2 * dz1);
        double B = (dx1 * dz2) - (dx2 * dz1);
        double C = (dx1 * dy2) - (dx2 * dy1);
        
        _a =  A;
        _b = -B;
        _c =  C;
        _d = A * p0[0] - B * p0[1] + C * p0[2];
        
        double norm = sqrt(_a*_a + _b*_b + _c*_c);
        _a /= norm;
        _b /= norm;
        _c /= norm;
        _d /= norm;
    }
    
    static Plane NotNormalizedPlane(Vec3d p0, Vec3d p1, Vec3d p2)
    {
        double dx1 = p1[0] - p0[0];
        double dx2 = p2[0] - p0[0];
        double dy1 = p1[1] - p0[1];
        double dy2 = p2[1] - p0[1];
        double dz1 = p1[2] - p0[2];
        double dz2 = p2[2] - p0[2];
        
        double A = (dy1 * dz2) - (dy2 * dz1);
        double B = (dx1 * dz2) - (dx2 * dz1);
        double C = (dx1 * dy2) - (dx2 * dy1);
        
        Plane alpha = Plane();
        
        alpha._a =  A;
        alpha._b = -B;
        alpha._c =  C;
        alpha._d = A * p0[0] - B * p0[1] + C * p0[2];
        
        return alpha;
    }
    
    
    Plane() // initialize a null plane
    {
        _a = 0.f;
        _b = 0.f;
        _c = 0.f;
        _d = 0.f;
    }
    
    inline double A() { return _a; }
    inline double B() { return _b; }
    inline double C() { return _c; }
    inline double D() { return _d; }
    
    inline Vec3d normal() { return Vec3d(_a, _b, _c); }
    
    inline bool isNull() { return (normal().length() == 0); }
    
    inline void invertOrientation()
    {
        _a *= -1;
        _b *= -1;
        _c *= -1;
        _d *= -1;
    }
    
    inline bool OnPlane( Vec3d point ) { return fabs(TestPoint(point)) < PLANE_EPSILON; }
    inline bool OnPositiveSide( Vec3d point )
    {
        double v = TestPoint(point);
        return (( fabs(v) > PLANE_EPSILON ) and ( v > 0.0 ));
    }
    inline bool OnNegativeSide( Vec3d point )
    {
        double v = TestPoint(point);
        return (( fabs(v) > PLANE_EPSILON ) and ( v < 0.0 ));
    }
    
    Vec3d ortho( Vec3d point_on_plane, Vec3d point )
    {
        Vec3d n     = normal();
        Vec3d v     = point - point_on_plane;
        double dist = dot(v,  n);
        n.normalize();
        Vec3d p     = point - n * dist;
        return p;
    }
    
private:
    double  _a, _b, _c, _d;
    double TestPoint( Vec3d p ) { return (p[0] * _a + p[1] * _b + p[2] * _c - _d); }
    
};

}}

#endif /* defined(__MeshEditE__Plane__) */
