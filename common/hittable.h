#ifndef HITTABLE_H
#define HITTABLE_H

#include "helpers.h"

class material;

struct hit_record
{
    point3 p;
    Vec3f normal;
    shared_ptr<material> mat_ptr;
    double t;
    bool front_face;

    inline void set_face_normal(const ray &r, const Vec3f &outward_normal)
    {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

#endif
