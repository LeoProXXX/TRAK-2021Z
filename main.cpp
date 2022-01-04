#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
//#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "common/external/tiny_obj_loader.h"

#include "common/helpers.h"

#include "common/camera.h"
#include "common/color.h"
#include "hittable_list.h"
#include "material.h"
#include "triangle.h"

#include <cstring>
#include <iostream>

color ray_color(const ray& r, const hittable& world, camera cam, int depth) {
    hit_record rec;

    std::map<std::string, vec3> light{
        {"position", vec3(5,0,-5)},
        {"ambient", vec3(0.3, 0.3, 0.3)},
        {"diffuse", vec3(1, 1, 1)},
        {"specular", vec3(1, 1, 1)}};
    

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0.001, infinity, rec))
    {
        color illumination = color(0, 0, 0);

        vec3 shiffted_point = rec.p + (1e-5 * rec.normal);
        vec3 intersection_to_light = unit_vector(light["position"] - shiffted_point);
        ray rr = ray(shiffted_point, intersection_to_light);

        hit_record rec_light;

        bool is_hitted = world.hit(rr, 0.001, infinity, rec_light);

        double min_distance = infinity;
        if (is_hitted)
        {
            min_distance = rec_light.t;
        }

        float intersection_to_light_distance = (light["position"] - rec.p).length();
        bool is_shadowed = min_distance < intersection_to_light_distance;

        if (is_shadowed)
        {
            return color(0, 0, 0);
        }

        // ambient
        illumination += rec.mat_ptr->ambient * light["ambient"];

        // diffuse
        illumination += rec.mat_ptr->diffuse * light["diffuse"] * dot(intersection_to_light, rec.normal);

        // specular
        vec3 intersection_to_camera = unit_vector(cam.origin - rec.p);
        vec3 H = unit_vector(intersection_to_light + intersection_to_camera);
        illumination += rec.mat_ptr->specular * light["specular"] * pow(dot(rec.normal, H), rec.mat_ptr->shininess);

        vec3 direction = reflect(r.direction(), rec.normal);
        ray new_ray = ray(shiffted_point, direction);

        return illumination + ray_color(new_ray, world, cam, depth-1);
    }

    return color(0, 0, 0);
}

// https://github.com/tinyobjloader/tinyobjloader
hittable_list loadOBJ(std::string input_file)
{
    hittable_list world;

    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./"; // Path to material

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
        std::cout << "TinyObjReader: " << reader.Warning();
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();
    auto &materials = reader.GetMaterials();

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
        {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + 0];
            tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
            tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
            tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

            vec3 v1(vx, vy, vz);

            idx = shapes[s].mesh.indices[index_offset + 1];
            vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
            vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
            vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

            vec3 v2(vx, vy, vz);

            idx = shapes[s].mesh.indices[index_offset + 2];
            vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
            vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
            vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

            vec3 v3(vx, vy, vz);

            tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
            tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
            tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];

            vec3 normal(nx, ny, nz);

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

            world.add(make_shared<triangle>(v1, v2, v3, normal, mat));

            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }
    }

    return world;
}

int main()
{
    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 1;
    const int max_depth = 50;

    // World - Read our .obj file
    hittable_list world = loadOBJ("cube.obj");

    // Camera
    point3 lookfrom(13,2,3);
    point3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

    // Render
    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            auto u = double(i) / (image_width-1);
            auto v = double(j) / (image_height-1);

            ray r = cam.get_ray(u, v);
            color pixel_color = ray_color(r, world, cam, max_depth);

            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
