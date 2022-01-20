#ifndef GRID_H
#define GRID_H

#include "../common/helpers.h"
#include "../common/mesh.h"
#include "acceleration_structure.h"

class Grid : public AccelerationStructure
{
    struct Cell
    {
        Cell() {}
        struct TriangleDesc
        {
            TriangleDesc(const Mesh *m, const uint32_t &t) : mesh(m), tri(t) {}
            const Mesh *mesh;
            uint32_t tri;
        };

        void insert(const Mesh* mesh, uint32_t t);
        bool intersect(const Vec3f &, const Vec3f &, const uint32_t &, float &, const Mesh *&, hit_record &rec) const;
        std::vector<TriangleDesc> triangles;
    };

public:
    Grid(std::vector<std::unique_ptr<const Mesh>> &m, uint32_t resolution_dim);
    ~Grid();
    bool intersect(const Vec3f &, const Vec3f &, const uint32_t &, float &, hit_record &) const;
    Cell **cells;
    BBox<> bbox;
    Vec3<uint32_t> resolution;
    Vec3f cellDimension;
};

#endif
