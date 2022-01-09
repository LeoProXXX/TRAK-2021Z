#ifndef CAMERA_H
#define CAMERA_H

#include "helpers.h"


class camera {
    public:
        camera() : camera(point3(0,0,-1), point3(0,0,0), Vec3f(0,1,0), 40, 1, 0, 10) {}

        camera(
            point3 lookfrom,
            point3 lookat,
            Vec3f   vup,
            double vfov, // vertical field-of-view in degrees
            double aspect_ratio,
            double aperture,
            double focus_dist,
            double _time0 = 0,
            double _time1 = 0
        ) {
            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta/2);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            w = lookfrom - lookat;
            normalize(w);
            // w = normalize(lookfrom - lookat);
            u = cross(vup, w);
            normalize(u);
            // u = normalize(cross(vup, w));
            v = cross(w, u);

            origin = lookfrom;
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;
            lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist*w;

            lens_radius = aperture / 2;
            time0 = _time0;
            time1 = _time1;
        }

        ray get_ray(double u, double v) const {
            return ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
        }

    public:
        point3 origin;
        point3 lower_left_corner;
        Vec3f horizontal;
        Vec3f vertical;
        Vec3f u, v, w;
        double lens_radius;
        double time0, time1;  // shutter open/close times
};

#endif
