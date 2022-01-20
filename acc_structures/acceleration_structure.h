#ifndef ACCELERATION_STRUCTURE_H
#define ACCELERATION_STRUCTURE_H

#include "../common/helpers.h"
#include "../common/mesh.h"

class AccelerationStructure
{
public:
    AccelerationStructure(std::vector<std::unique_ptr<const Mesh>> &m) : meshes(std::move(m)) {}
    virtual ~AccelerationStructure() {}
    virtual bool intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit, hit_record& rec) const;

protected:
    const std::vector<std::unique_ptr<const Mesh>> meshes;
};

#endif
