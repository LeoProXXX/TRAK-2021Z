#ifndef BBOX_H
#define BBOX_H

#include "helpers.h"

template <typename T = float>
class BBox
{
public:
    BBox() {}
    BBox(Vec3<T> min_, Vec3<T> max_)
    {
        bounds[0] = min_;
        bounds[1] = max_;
    }
    BBox &extendBy(const Vec3<T> &p)
    {
        if (p.x < bounds[0].x)
            bounds[0].x = p.x;
        if (p.y < bounds[0].y)
            bounds[0].y = p.y;
        if (p.z < bounds[0].z)
            bounds[0].z = p.z;
        if (p.x > bounds[1].x)
            bounds[1].x = p.x;
        if (p.y > bounds[1].y)
            bounds[1].y = p.y;
        if (p.z > bounds[1].z)
            bounds[1].z = p.z;

        return *this;
    }
    /*inline */ Vec3<T> centroid() const { return (bounds[0] + bounds[1]) * 0.5; }
    Vec3<T> &operator[](bool i) { return bounds[i]; }
    const Vec3<T> operator[](bool i) const { return bounds[i]; }
    bool intersect(const Vec3<T> &, const Vec3<T> &, const Vec3b &, float &) const;
    Vec3<T> bounds[2] = {kInfinity, -kInfinity};
};

template <typename T>
bool BBox<T>::intersect(const Vec3<T> &orig, const Vec3<T> &invDir, const Vec3b &sign, float &tHit) const
{
    numRayBBoxTests++;
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bounds[sign[0]].x - orig.x) * invDir.x;
    tmax = (bounds[1 - sign[0]].x - orig.x) * invDir.x;
    tymin = (bounds[sign[1]].y - orig.y) * invDir.y;
    tymax = (bounds[1 - sign[1]].y - orig.y) * invDir.y;

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    tzmin = (bounds[sign[2]].z - orig.z) * invDir.z;
    tzmax = (bounds[1 - sign[2]].z - orig.z) * invDir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    tHit = tmin;

    return true;
}

#endif
