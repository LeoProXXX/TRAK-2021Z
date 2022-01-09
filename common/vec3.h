#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

// using std::sqrt;
// using std::fabs;

// class vec3 {
//     public:
//         vec3() : e{0,0,0} {}
//         vec3(double ee) : e{ee, ee, ee} {}
//         vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

//         double x() const { return e[0]; }
//         double y() const { return e[1]; }
//         double z() const { return e[2]; }

//         vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
//         double operator[](int i) const { return e[i]; }
//         double& operator[](int i) { return e[i]; }


//         vec3 operator / (const vec3& v) const { return vec3(this->x() / v.x(), this->y() / v.y(), this->z() / v.z()); }
//         // friend vec3 operator / (const T r, const Vec3& v)
//         // { return Vec3(r / v.x, r / v.y, r / v.z); }

//         vec3& operator+=(const vec3 &v) {
//             e[0] += v.e[0];
//             e[1] += v.e[1];
//             e[2] += v.e[2];
//             return *this;
//         }

//         vec3& operator*=(const double t) {
//             e[0] *= t;
//             e[1] *= t;
//             e[2] *= t;
//             return *this;
//         }

//         vec3& operator/=(const double t) {
//             return *this *= 1/t;
//         }

//         double length() const {
//             return sqrt(length_squared());
//         }

//         double length_squared() const {
//             return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
//         }

//         bool near_zero() const {
//             // Return true if the vector is close to zero in all dimensions.
//             const auto s = 1e-8;
//             return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
//         }

//         inline static vec3 random() {
//             return vec3(random_double(), random_double(), random_double());
//         }

//         inline static vec3 random(double min, double max) {
//             return vec3(random_double(min,max), random_double(min,max), random_double(min,max));
//         }

//     public:
//         double e[3];
// };


// Type aliases for vec3
// using point3 = vec3;   // 3D point
// using color = vec3;    // RGB color


// // vec3 Utility Functions

// inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
//     return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
// }

// inline vec3 operator+(const vec3 &u, const vec3 &v) {
//     return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
// }

// inline vec3 operator-(const vec3 &u, const vec3 &v) {
//     return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
// }

// inline vec3 operator*(const vec3 &u, const vec3 &v) {
//     return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
// }

// inline vec3 operator*(double t, const vec3 &v) {
//     return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
// }

// inline vec3 operator*(const vec3 &v, double t) {
//     return t * v;
// }

// inline vec3 operator/(vec3 v, double t) {
//     return (1/t) * v;
// }

// inline double dot(const vec3 &u, const vec3 &v) {
//     return u.e[0] * v.e[0]
//          + u.e[1] * v.e[1]
//          + u.e[2] * v.e[2];
// }

// inline vec3 cross(const vec3 &u, const vec3 &v) {
//     return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
//                 u.e[2] * v.e[0] - u.e[0] * v.e[2],
//                 u.e[0] * v.e[1] - u.e[1] * v.e[0]);
// }

// inline vec3 unit_vector(vec3 v) {
//     return v / v.length();
// }

// inline vec3 random_in_unit_disk() {
//     while (true) {
//         auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
//         if (p.length_squared() >= 1) continue;
//         return p;
//     }
// }

// inline vec3 random_in_unit_sphere() {
//     while (true) {
//         auto p = vec3::random(-1,1);
//         if (p.length_squared() >= 1) continue;
//         return p;
//     }
// }

// inline vec3 random_unit_vector() {
//     return unit_vector(random_in_unit_sphere());
// }

// inline vec3 random_in_hemisphere(const vec3& normal) {
//     vec3 in_unit_sphere = random_in_unit_sphere();
//     if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
//         return in_unit_sphere;
//     else
//         return -in_unit_sphere;
// }

// inline vec3 reflect(const vec3& v, const vec3& n) {
//     return v - 2*dot(v,n)*n;
// }

// inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
//     auto cos_theta = fmin(dot(-uv, n), 1.0);
//     vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
//     vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
//     return r_out_perp + r_out_parallel;
// }




template<typename T>
class Vec3
{
public:
    Vec3() : x(0), y(0), z(0) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3 operator * (const T& r) const { return Vec3(x * r, y * r, z * r); }
    Vec3 operator + (const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator - (const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }

    Vec3 &operator+=(const Vec3 &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    template<typename U>
    Vec3 operator / (const Vec3<U>& v) const { return Vec3(x / v.x, y / v.y, z / v.z); }
    friend Vec3 operator / (const T r, const Vec3& v)
    { return Vec3(r / v.x, r / v.y, r / v.z); }
    Vec3 operator - () const { return Vec3(-x, -y, -z); }
    const T& operator [] (size_t i) const { return (&x)[i]; }
    T& operator [] (size_t i) { return (&x)[i]; }
    T length2() const{ return x * x + y * y + z * z; }
    T length() const{ return sqrt(x * x + y * y + z * z); }
    friend Vec3 operator * (const float&r, const Vec3& v)
    { return Vec3(v.x * r, v.y * r, v.z * r); }
    friend std::ostream& operator << (std::ostream& os, const Vec3<T>& v)
    { os << v.x << " " << v.y << " " << v.z << std::endl; return os; }
	T x, y, z;
};

template<typename T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
    return Vec3<T>(a.y * b.z - a.z * b.y,
                   a.z * b.x - a.x * b.z,
                   a.x * b.y - a.y * b.x);
}

template<typename T>
T dot(const Vec3<T>& va, const Vec3<T>& vb)
{ return va.x * vb.x + va.y * vb.y + va.z * vb.z; }

template<typename T>
void normalize(Vec3<T>& vec)
{
    T len2 = vec.length2();
    // if (len2 > 0) {
        T invLen = 1 / sqrt(len2);
        vec.x *= invLen, vec.y *= invLen, vec.z *= invLen;
    // }
}

template<typename T>
class Matrix44
{
public:
    Matrix44() { /* ... define identity matrix ... */ }
    Matrix44(T m00, T m01, T m02, T m03,
             T m10, T m11, T m12, T m13,
             T m20, T m21, T m22, T m23,
             T m30, T m31, T m32, T m33)
    {
        m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
        m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
        m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
        m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
    }
    Matrix44 inverse() const { Matrix44 matInv = *this; return matInv; }
    T* operator [] (size_t i) { return &m[i][0]; }
    const T* operator [] (size_t i) const { return &m[i][0]; }
    T m[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
};

template<typename T>
void matVecMult(const Matrix44<T>& m, Vec3<T>& v)
{
    Vec3<T> vt;
    vt.x = v.x * m[0][0] + v.y * m[1][0] + v.z * m[2][0],
    vt.y = v.x * m[0][1] + v.y * m[1][1] + v.z * m[2][1],
    vt.z = v.x * m[0][2] + v.y * m[1][2] + v.z * m[2][2];

    v = vt;
}

template<typename T>
void matPointMult(const Matrix44<T>& m, Vec3<T>& p)
{
    Vec3<T> pt;
    pt.x = p.x * m[0][0] + p.y * m[1][0] + p.z * m[2][0] + m[3][0];
    pt.y = p.x * m[0][1] + p.y * m[1][1] + p.z * m[2][1] + m[3][1];
    pt.z = p.x * m[0][2] + p.y * m[1][2] + p.z * m[2][2] + m[3][2];
    T w  = p.x * m[0][3] + p.y * m[1][3] + p.z * m[2][3] + m[3][3];
    if (w != 1) {
        pt.x /= w, pt.y /= w, pt.z /= w;
    }

    p = pt;
}

using Vec3f = Vec3<float>;
using Vec3b = Vec3<bool>;
using Vec3i = Vec3<int32_t>;
using Vec3ui = Vec3<uint32_t>;
using Matrix44f = Matrix44<float>;

using point3 = Vec3f;   // 3D point
using color = Vec3f;    // RGB color

inline Vec3f operator/(Vec3f v, float t) {
    return (1/t) * v;
}

inline Vec3f reflect(const Vec3f& v, const Vec3f& n) {
    return v - 2*dot(v,n)*n;
}

inline Vec3f operator*(const Vec3f &u, const Vec3f &v) {
    return Vec3f(u.x * v.x, u.y * v.y, u.z * v.z);
}


#endif
