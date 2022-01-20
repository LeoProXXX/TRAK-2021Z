#ifndef MESH_H
#define MESH_H

#include "helpers.h"
#include "bbox.h"
#include "hittable.h"
#include "material.h"

class GeomPrimitive
{
public:
    GeomPrimitive(Matrix44f objectToWorld_) : objectToWorld(objectToWorld_), worldToObject(objectToWorld.inverse())
    {
    }
    virtual ~GeomPrimitive() {}
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, hit_record &) const = 0;
    Matrix44f objectToWorld;
    Matrix44f worldToObject;
    BBox<> bbox;
    uint32_t test;
};

class Mesh : public GeomPrimitive
{
public:
    Mesh(
        uint32_t numPolygons,
        const std::vector<uint32_t> &polygonNumVertsArray,
        const std::vector<uint32_t> &polygonIndicesInVertexPool,
        std::vector<Vec3f> vertexPool_,
        std::vector<Vec3f> normals_,
        std::vector<shared_ptr<material>> matrials_,
        Matrix44f objectToWorld_ = Matrix44f()) : GeomPrimitive(objectToWorld_),
                                                  vertexPool(vertexPool_),
                                                  matrials(matrials_),
                                                  normals(normals_)
    {
        for (uint32_t i = 0; i < vertexPool.size(); ++i)
        {
            matPointMult(objectToWorld, vertexPool[i]);
            bbox.extendBy(vertexPool[i]);
        }
        // compute total number of triangles
        for (uint32_t i = 0; i < numPolygons; ++i)
        {
            assert(polygonNumVertsArray[i] >= 3);
            numTriangles += polygonNumVertsArray[i] - 2;
        }

        triangleIndicesInVertexPool.reserve(numTriangles * 3);
        mailbox.reserve(numTriangles);
        // for each face
        for (uint32_t i = 0, offset = 0, currTriangleIndex = 0; i < numPolygons; ++i)
        {
            // for each triangle in the face
            for (uint32_t j = 0; j < polygonNumVertsArray[i] - 2; ++j)
            {
                triangleIndicesInVertexPool[currTriangleIndex] = polygonIndicesInVertexPool[offset];
                triangleIndicesInVertexPool[currTriangleIndex + 1] = polygonIndicesInVertexPool[offset + j + 1];
                triangleIndicesInVertexPool[currTriangleIndex + 2] = polygonIndicesInVertexPool[offset + j + 2];
                currTriangleIndex += 3;
            }
            offset += polygonNumVertsArray[i];
        }
    }
    bool intersect(const Vec3f &rayOrig, const Vec3f &rayDir, float &tNear, hit_record &rec) const;

    uint32_t numTriangles = {0};
    std::vector<uint32_t> triangleIndicesInVertexPool;
    std::vector<Vec3f> vertexPool;

    std::vector<Vec3f> normals;
    std::vector<shared_ptr<material>> matrials;

    mutable std::vector<uint32_t> mailbox;
};

// use Moller-Trumbor method
inline bool rayTriangleIntersect(
    const Vec3f &orig, const Vec3f &dir,
    const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
    float &t, float &u, float &v)
{
    numRayTriangleTests++;
    Vec3f v0v1 = v1 - v0;
    Vec3f v0v2 = v2 - v0;
    Vec3f pvec = cross(dir, v0v2);
    float det = dot(v0v1, pvec);

    // ray and triangle are parallel if det is close to 0
    if (det < kEpsilon)
        return false;
    // if (fabs(det) < kEpsilon) return false;

    float invDet = 1 / det;

    Vec3f tvec = orig - v0;
    u = dot(tvec, pvec) * invDet;

    if (u < 0 || u > 1)
        return false;

    Vec3f qvec = cross(tvec, v0v1);
    v = dot(dir, qvec) * invDet;

    if (v < 0 || u + v > 1)
        return false;

    t = dot(v0v2, qvec) * invDet;

    if (t < 0)
        return false;

    numRayTriangleIntersections++;

    return true;
}

inline bool Mesh::intersect(const Vec3f &rayOrig, const Vec3f &rayDir, float &tNear, hit_record &rec) const
{
    float t, u, v;
    uint32_t intersectedTriIndex;
    bool intersected = false;

    for (uint32_t i = 0; i < numTriangles; ++i)
    {
        if (rayTriangleIntersect(rayOrig, rayDir,
                                 vertexPool[triangleIndicesInVertexPool[i * 3]],
                                 vertexPool[triangleIndicesInVertexPool[i * 3 + 1]],
                                 vertexPool[triangleIndicesInVertexPool[i * 3 + 2]], t, u, v) &&
            t < tNear)
        {
            tNear = t;
            intersectedTriIndex = i;
            intersected = true;

            ray r = ray(rayOrig, rayDir);
            rec.t = t;
            rec.p = r.at(rec.t);
            rec.set_face_normal(r, normals[i]);
            rec.mat_ptr = matrials[i];
        }
    }

    return intersected;
}

#endif
