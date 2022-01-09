#ifndef RAY_H
#define RAY_H

#include "vec3.h"


class ray {
    public:
        ray() {}
        ray(const point3& origin, const Vec3f& direction)
            : orig(origin), dir(direction), tm(0)
        {}

        ray(const point3& origin, const Vec3f& direction, double time)
            : orig(origin), dir(direction), tm(time)
        {}

        point3 origin() const  { return orig; }
        Vec3f direction() const { return dir; }
        double time() const    { return tm; }

        point3 at(double t) const {
            return orig + t*dir;
        }

    public:
        point3 orig;
        Vec3f dir;
        double tm;
};

#endif
