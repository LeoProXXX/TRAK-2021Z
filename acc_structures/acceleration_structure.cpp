#include "acceleration_structure.h"


bool AccelerationStructure::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit, hit_record& rec) const
{
    const Mesh* intersectedMesh = nullptr;
    float t = kInfinity;
    for (const auto& mesh : meshes)
    {
        if (mesh->intersect(orig, dir, t, rec) && t < tHit)
        {
            intersectedMesh = mesh.get();
            tHit = t;
        }
    }

    return (intersectedMesh != nullptr);
}
