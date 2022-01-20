#ifndef MESHES_DATA_PROVIDER_H
#define MESHES_DATA_PROVIDER_H

#include "helpers.h"
#include "mesh.h"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
//#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "external/tiny_obj_loader.h"

// https://github.com/tinyobjloader/tinyobjloader
std::vector<std::unique_ptr<const Mesh>> loadOBJMeshes(std::string input_file)
{
    std::vector<std::unique_ptr<const Mesh>> meshes;

    tinyobj::ObjReaderConfig reader_config;
    // reader_config.mtl_search_path = "./"; // Path to material

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(input_file, reader_config))
    {
        if (!reader.Error().empty())
        {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty())
    {
        std::cerr << "TinyObjReader: " << reader.Warning();
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();
    auto &materials = reader.GetMaterials();

    uint32_t numPolygons = 0;

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
        numPolygons += shapes[s].mesh.num_face_vertices.size();
    }

    std::vector<uint32_t> polyNumVertsArray(numPolygons, 3);
    std::vector<uint32_t> polyIndicesInVertPool(numPolygons * 3);

    std::vector<Vec3f> normals(numPolygons);
    std::vector<shared_ptr<material>> materialss(numPolygons);

    std::vector<Vec3f> vertPool((int)(attrib.vertices.size() / 3));

    for (size_t i = 0; i < vertPool.size(); i++)
    {
        vertPool[i] = Vec3f(attrib.vertices[3 * size_t(i) + 0], attrib.vertices[3 * size_t(i) + 1], attrib.vertices[3 * size_t(i) + 2]);
    }

    int offset = 0;
    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + 0];
            polyIndicesInVertPool[offset] = idx.vertex_index;

            idx = shapes[s].mesh.indices[index_offset + 1];
            polyIndicesInVertPool[offset + 1] = idx.vertex_index;

            idx = shapes[s].mesh.indices[index_offset + 2];
            polyIndicesInVertPool[offset + 2] = idx.vertex_index;

            tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
            tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
            tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
            Vec3f normal(nx, ny, nz);
            normals[offset / 3] = normal;

            auto mat = make_shared<material>(
                color(
                    materials[shapes[s].mesh.material_ids[f]].ambient[0],
                    materials[shapes[s].mesh.material_ids[f]].ambient[1],
                    materials[shapes[s].mesh.material_ids[f]].ambient[2]),
                color(
                    materials[shapes[s].mesh.material_ids[f]].diffuse[0],
                    materials[shapes[s].mesh.material_ids[f]].diffuse[1],
                    materials[shapes[s].mesh.material_ids[f]].diffuse[2]),
                color(
                    materials[shapes[s].mesh.material_ids[f]].specular[0],
                    materials[shapes[s].mesh.material_ids[f]].specular[1],
                    materials[shapes[s].mesh.material_ids[f]].specular[2]),
                materials[shapes[s].mesh.material_ids[f]].shininess);

            materialss[offset / 3] = mat;

            offset += 3;

            index_offset += fv;
        }
    }

    meshes.emplace_back(new Mesh(numPolygons, polyNumVertsArray, polyIndicesInVertPool, vertPool, normals, materialss));

    return meshes;
}

#endif
