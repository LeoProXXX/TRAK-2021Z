#include "grid.h"
#include <cstring>


Grid::Grid(std::vector<std::unique_ptr<const Mesh>>& m) : AccelerationStructure(m)
{
    uint32_t totalNumTriangles = 0;
    for (const auto& m : meshes)
    {
        bbox.extendBy(m->bbox[0]);
        bbox.extendBy(m->bbox[1]);
        totalNumTriangles += m->numTriangles;
    }
    // Create the grid
    Vec3f size = bbox[1] - bbox[0];
    float cubeRoot = powf(totalNumTriangles / (size.x * size.y * size.z), 1. / 3.f);
    for (uint8_t i = 0; i < 3; ++i)
    {
        resolution[i] = std::floor(size[i] * cubeRoot);
        if (resolution[i] < 1)
            resolution[i] = 1;
        if (resolution[i] > 128)
            resolution[i] = 128;
    }
    cellDimension = size / resolution;

    uint32_t numCells = resolution.x * resolution.y * resolution.z;
    cells = new Grid::Cell * [numCells];
    memset(cells, 0x0, sizeof(Grid::Cell*) * numCells);

    for (const auto& m : meshes)
    {
        for (uint32_t i = 0, off = 0; i < m->numTriangles; ++i, off += 3)
        {
            Vec3f min(kInfinity), max(-kInfinity);
            const Vec3f& v0 = m->vertexPool[m->triangleIndicesInVertexPool[off]];
            const Vec3f& v1 = m->vertexPool[m->triangleIndicesInVertexPool[off + 1]];
            const Vec3f& v2 = m->vertexPool[m->triangleIndicesInVertexPool[off + 2]];

            for (uint8_t j = 0; j < 3; ++j)
            {
                if (v0[j] < min[j])
                    min[j] = v0[j];
                if (v1[j] < min[j])
                    min[j] = v1[j];
                if (v2[j] < min[j])
                    min[j] = v2[j];
                if (v0[j] > max[j])
                    max[j] = v0[j];
                if (v1[j] > max[j])
                    max[j] = v1[j];
                if (v2[j] > max[j])
                    max[j] = v2[j];
            }
            // Convert to cell coordinates
            min = (min - bbox[0]) / cellDimension;
            max = (max - bbox[0]) / cellDimension;
            uint32_t zmin = clamp<int32_t>(std::floor(min[2]), 0, resolution[2] - 1);
            uint32_t zmax = clamp<int32_t>(std::floor(max[2]), 0, resolution[2] - 1);
            uint32_t ymin = clamp<int32_t>(std::floor(min[1]), 0, resolution[1] - 1);
            uint32_t ymax = clamp<int32_t>(std::floor(max[1]), 0, resolution[1] - 1);
            uint32_t xmin = clamp<int32_t>(std::floor(min[0]), 0, resolution[0] - 1);
            uint32_t xmax = clamp<int32_t>(std::floor(max[0]), 0, resolution[0] - 1);
            // Loop over the cells the triangle overlaps and insert
            for (uint32_t z = zmin; z <= zmax; ++z)
            {
                for (uint32_t y = ymin; y <= ymax; ++y)
                {
                    for (uint32_t x = xmin; x <= xmax; ++x)
                    {
                        uint32_t index = z * resolution[0] * resolution[1] + y * resolution[0] + x;
                        if (cells[index] == NULL)
                            cells[index] = new Grid::Cell;
                        cells[index]->insert(m.get(), i);
                    }
                }
            }
        }
    }
}

Grid::~Grid()
{
    for (uint32_t i = 0; i < resolution[0] * resolution[1] * resolution[2]; ++i)
        if (cells[i] != NULL)
            delete cells[i];
    delete[] cells;
}

void Grid::Cell::insert(const Mesh* mesh, uint32_t t)
{
    triangles.push_back(Grid::Cell::TriangleDesc(mesh, t));
}

bool Grid::Cell::intersect(
    const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId,
    float& tHit, const Mesh*& intersectedMesh, hit_record& rec) const
{
    float uhit, vhit;
    for (uint32_t i = 0; i < triangles.size(); ++i)
    {
        if (rayId != triangles[i].mesh->mailbox[triangles[i].tri])
        {
            triangles[i].mesh->mailbox[triangles[i].tri] = rayId;
            const Mesh* mesh = triangles[i].mesh;
            uint32_t j = triangles[i].tri * 3;
            const Vec3f& v0 = mesh->vertexPool[mesh->triangleIndicesInVertexPool[j]];
            const Vec3f& v1 = mesh->vertexPool[mesh->triangleIndicesInVertexPool[j + 1]];
            const Vec3f& v2 = mesh->vertexPool[mesh->triangleIndicesInVertexPool[j + 2]];
            float t, u, v;

            if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v))
            {
                if (t < tHit)
                {
                    tHit = t;
                    uhit = u;
                    vhit = v;
                    intersectedMesh = triangles[i].mesh;

                    ray r = ray(orig, dir);
                    rec.t = t;
                    rec.p = r.at(rec.t);
                    rec.set_face_normal(r, triangles[i].mesh->normals[j / 3]);
                    rec.mat_ptr = triangles[i].mesh->matrials[j / 3];
                }
            }
        }
    }

    return (intersectedMesh != nullptr);
}


bool Grid::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit, hit_record& rec) const
{
    const Vec3f invDir = 1 / dir;
    const Vec3b sign(dir.x < 0, dir.y < 0, dir.z < 0);
    float tHitBox;
    if (!bbox.intersect(orig, invDir, sign, tHitBox))
        return false;

    // initialization step
    Vec3i exit, step, cell;
    Vec3f deltaT, nextCrossingT;
    for (uint8_t i = 0; i < 3; ++i)
    {
        float rayOrigCell = ((orig[i] + dir[i] * tHitBox) - bbox[0][i]);
        cell[i] = clamp<int32_t>(std::floor(rayOrigCell / cellDimension[i]), 0, resolution[i] - 1);
        if (dir[i] < 0)
        {
            deltaT[i] = -cellDimension[i] * invDir[i];
            nextCrossingT[i] = tHitBox + (cell[i] * cellDimension[i] - rayOrigCell) * invDir[i];
            exit[i] = -1;
            step[i] = -1;
        }
        else
        {
            deltaT[i] = cellDimension[i] * invDir[i];
            nextCrossingT[i] = tHitBox + ((cell[i] + 1) * cellDimension[i] - rayOrigCell) * invDir[i];
            exit[i] = resolution[i];
            step[i] = 1;
        }
    }

    const Mesh* intersectedMesh = nullptr;
    while (1)
    {
        uint32_t o = cell[2] * resolution[0] * resolution[1] + cell[1] * resolution[0] + cell[0];
        if (cells[o] != nullptr)
        {
            cells[o]->intersect(orig, dir, rayId, tHit, intersectedMesh, rec);
        }
        uint8_t k =
            ((nextCrossingT[0] < nextCrossingT[1]) << 2) +
            ((nextCrossingT[0] < nextCrossingT[2]) << 1) +
            ((nextCrossingT[1] < nextCrossingT[2]));
        static const uint8_t map[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
        uint8_t axis = map[k];

        if (tHit < nextCrossingT[axis])
            break;
        cell[axis] += step[axis];
        if (cell[axis] == exit[axis])
            break;
        nextCrossingT[axis] += deltaT[axis];
    }

    return (intersectedMesh != nullptr);
}
