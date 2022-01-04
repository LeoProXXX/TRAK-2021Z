#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "common/helpers.h"

#include "hittable.h"

constexpr double kEpsilon = 1e-8;

class triangle : public hittable
{
public:
    triangle() {}

    triangle(vec3 vv1, vec3 vv2, vec3 vv3, vec3 normall, shared_ptr<material> m) : v1(vv1), v2(vv2), v3(vv3), normal(normall), mat_ptr(m){};

    virtual bool hit(
        const ray &r, double t_min, double t_max, hit_record &rec) const override;

public:
    vec3 v1;
    vec3 v2;
    vec3 v3;
    vec3 normal;
    shared_ptr<material> mat_ptr;
};

bool triangle::hit(const ray &r, double t_min, double t_max, hit_record &rec) const
{
    // const float EPSILON = 0.0000001;
    // vec3 vertex0 = v1;
    // vec3 vertex1 = v2;  
    // vec3 vertex2 = v3;
    // vec3 edge1, edge2, h, s, q;
    // float a,f,u,v;
    // edge1 = vertex1 - vertex0;
    // edge2 = vertex2 - vertex0;
    // h = cross(r.direction(), edge2);
    // a = dot(edge1, h);
    // if (a > -EPSILON && a < EPSILON)
    //     return false;    // This ray is parallel to this triangle.
    // f = 1.0/a;
    // s = r.origin() - vertex0;
    // u = f * dot(s, h);
    // if (u < 0.0 || u > 1.0)
    //     return false;
    // q = cross(s, edge1);
    // v = f * dot(r.direction(), q);
    // if (v < 0.0 || u + v > 1.0)
    //     return false;
    // // At this stage we can compute t to find out where the intersection point is on the line.
    // float t = f * dot(edge2, q);
    // if (t > EPSILON) // ray intersection
    // {
    //     rec.t = t;
    //     rec.p = r.at(rec.t);
    //     rec.set_face_normal(r, normal);
    //     rec.mat_ptr = mat_ptr;
    //     return true;
    // }
    // else // This means that there is a line intersection but not a ray intersection.
    //     return false;



    // https://www.scratchapixel.com/code.php?id=9&origin=/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle

#define MOLLER_TRUMBORE
#define CULLING

#ifdef MOLLER_TRUMBORE 
    vec3 v0v1 = v2 - v1; 
    vec3 v0v2 = v3 - v1; 
    vec3 pvec = cross(r.direction(), v0v2); 
    double det = dot(v0v1, pvec); 

#ifdef CULLING 
    if (det < kEpsilon) return false;  
#else 
    // ray and triangle are parallel if det is close to 0
    if (fabs(det) < kEpsilon) return false;
#endif 
    double invDet = 1 / det;

    vec3 tvec = r.origin() - v1;
    double u = dot(tvec, pvec) * invDet;
    if (u < 0 || u > 1) return false;

    vec3 qvec = cross(tvec, v0v1);
    double v = dot(r.direction(), qvec) * invDet; 
    if (v < 0 || u + v > 1) return false; 

    double t = dot(v0v2, qvec) * invDet;

    if (t < 0 || t > t_max) return false;

    rec.t = t;
    rec.p = r.at(rec.t);
    rec.set_face_normal(r, normal);
    rec.mat_ptr = mat_ptr;

    return true;

#else
    // compute plane's normal
    vec3 v0v1 = v2 - v1; 
    vec3 v0v2 = v3 - v1;
    // no need to normalize
    vec3 N = cross(v0v1, v0v2); 
    float denom = dot(N, N);
 
    // Step 1: finding P
 
    // check if ray and plane are parallel ?
    double NdotRayDirection = dot(N, r.direction()); 
    if (fabs(NdotRayDirection) < kEpsilon) // almost 0 
        return false; // they are parallel so they don't intersect ! 
 
    // compute d parameter using equation 2
    double d = dot(N, v1);
 
    // compute t (equation 3)
    double t = (dot(N, r.origin()) + d) / NdotRayDirection; 
    // check if the triangle is in behind the ray
    if (t < 0) return false; // the triangle is behind
 
    // compute the intersection point using equation 1
    vec3 P = r.origin() + t * r.direction(); 
 
    // Step 2: inside-outside test
    vec3 C; // vector perpendicular to triangle's plane 
 
    // edge 0
    vec3 edge0 = v2 - v1; 
    vec3 vp0 = P - v1; 
    C = cross(edge0, vp0); 
    if (dot(N, C) < 0) return false; // P is on the right side 
 
    // edge 1
    float u;
    vec3 edge1 = v2 - v1; 
    vec3 vp1 = P - v1; 
    C = cross(edge1, vp1); 
    if ((u = dot(N, C)) < 0)  return false; // P is on the right side 
 
    // edge 2
    float v;
    vec3 edge2 = v1 - v3; 
    vec3 vp2 = P - v3; 
    C = cross(edge2, vp2); 
    if ((v = dot(N, C)) < 0) return false; // P is on the right side; 
 
    u /= denom;
    v /= denom;

    rec.t = t;
    rec.p = r.at(rec.t);
    rec.set_face_normal(r, normal);
    rec.mat_ptr = mat_ptr;
 
    return true; // this ray hits the triangle
#endif
}

#endif
