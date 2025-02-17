#ifndef HRLPERS_H
#define HRLPERS_H

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>
#include <cassert>
#include <queue>

#include <atomic>
inline std::atomic<uint32_t> numPrimaryRays(0);
inline std::atomic<uint32_t> numRayTriangleTests(0);
inline std::atomic<uint32_t> numRayTriangleIntersections(0);
inline std::atomic<uint32_t> numRayBBoxTests(0);
inline std::atomic<uint32_t> numRayBoundingVolumeTests(0);

// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;

// Constants
const float kInfinity = std::numeric_limits<float>::max();
constexpr double kEpsilon = 1e-8;
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}

inline int random_int(int min, int max) {
    // Returns a random integer in [min,max].
    return static_cast<int>(random_double(min, max+1));
}

template<typename T> inline T clamp(const T &v, const T &lo, const T &hi)
{ return std::max(lo, std::min(v, hi)); }

// Common Headers

#include "ray.h"
#include "vec3.h"


#endif
