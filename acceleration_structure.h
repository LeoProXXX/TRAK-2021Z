#ifndef ACCELERATION_STRUCTURE_H
#define ACCELERATION_STRUCTURE_H

#include "common/helpers.h"
#include "common/mesh.h"

// [comment]
// The most basic acceleration class (the parent class of all the other acceleration structures)
// could have a *pure* virtual intersect() method but instead we decided in this implementation
// to have it supporting the basic ray-mesh intersection routine.
// [/comment]
class AccelerationStructure
{
public:
    // [comment]
    // We transfer owner ship of the mesh list to the acceleration structure. This makes
    // more sense from a functional/structure stand point because the objects/meshes themselves
    // should be destroyed/deleted when the acceleration structure is being deleted
    // Ideally this means the render function() itself should be bounded (in terms of lifespan)
    // to the lifespan of the acceleration structure (aka we should wrap the accel struc instance
    // and the render method() within the same object, so that when this object is deleted,
    // the render function can't be called anymore.
    // [/comment]
    AccelerationStructure(std::vector<std::unique_ptr<const Mesh>> &m) : meshes(std::move(m)) {}
    virtual ~AccelerationStructure() {}
    virtual bool intersect(const Vec3f &orig, const Vec3f &dir, const uint32_t &rayId, float &tHit, hit_record &rec) const
    {
        // [comment]
        // Because we don't want to change the content of the mesh itself, just get a point to it so
        // it's safer to make it const (which doesn't mean we can't change its assignment just that
        // we can't do something like intersectedMesh->color = blue. You would get something like:
        // "read-only variable is not assignable" error message at compile time)
        // [/comment]
        const Mesh *intersectedMesh = nullptr;
        float t = kInfinity;
        for (const auto &mesh : meshes)
        {
            if (mesh->intersect(orig, dir, t, rec) && t < tHit)
            {
                intersectedMesh = mesh.get();
                tHit = t;
            }
        }

        return (intersectedMesh != nullptr);
    }

protected:
    const std::vector<std::unique_ptr<const Mesh>> meshes;
};

#endif
