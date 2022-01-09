#ifndef ACCELERATION_STRUCTURE_H
#define ACCELERATION_STRUCTURE_H

#include "common/helpers.h"
#include "common/mesh.h"

class AccelerationStructure
{
public:
    AccelerationStructure(std::vector<std::unique_ptr<const Mesh>> &m) : meshes(std::move(m)) {}
    virtual ~AccelerationStructure() {}
    virtual bool intersect(const Vec3f &orig, const Vec3f &dir, const uint32_t &rayId, float &tHit, hit_record &rec) const
    {
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
