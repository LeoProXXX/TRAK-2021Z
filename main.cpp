#include "common/helpers.h"
#include "common/camera.h"
#include "common/color.h"
#include "common/material.h"
#include "common/mesh.h"
#include "meshes_data_provider.h"
#include "acceleration_structure.h"
#include "grid.h"

#include <cstring>
#include <iostream>
#include <time.h>

long long int rayIdd = 1;

std::map<std::string, Vec3f> light{
    {"position", Vec3f(5, 0, -5)},
    {"ambient", Vec3f(0.3, 0.3, 0.3)},
    {"diffuse", Vec3f(1, 1, 1)},
    {"specular", Vec3f(1, 1, 1)}};

void makeScene(std::vector<std::unique_ptr<const Mesh>> &meshes, std::string model_file_name)
{
    meshes = std::move(loadOBJMeshes(model_file_name));
}

color ray_color(const ray &r, const std::unique_ptr<AccelerationStructure> &accel, camera cam, int depth)
{
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    float tHit = kInfinity;

    if (accel->intersect(r.origin(), r.direction(), rayIdd++, tHit, rec))
    {   
        color illumination = color(0, 0, 0);

        Vec3f shiffted_point = rec.p + (1e-5 * rec.normal);
        Vec3f intersection_to_light = light["position"] - shiffted_point;
        normalize(intersection_to_light);

        ray rr = ray(shiffted_point, intersection_to_light);

        hit_record rec_light;

        float tHit2 = kInfinity;
        bool is_hitted = accel->intersect(rr.origin(), rr.direction(), rayIdd++, tHit2, rec_light);

        float min_distance = kInfinity;
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
        Vec3f intersection_to_camera = cam.origin - rec.p;
        normalize(intersection_to_camera);
        Vec3f H = intersection_to_light + intersection_to_camera;
        normalize(H);
        illumination += rec.mat_ptr->specular * light["specular"] * pow(dot(rec.normal, H), rec.mat_ptr->shininess);

        Vec3f direction = reflect(r.direction(), rec.normal);
        ray new_ray = ray(shiffted_point, direction);

        return illumination + ray_color(new_ray, accel, cam, depth - 1);
    }

    return color(0, 0, 0);
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cout << "Please provide 2 arguments: 1st - input .obj file; 2nd - output .ppm file; 3rd - acc data structure: [none, grid, bvh]"
                  << std::endl;
        return 1;
    }

    time_t timer;
    float time_diff;

    std::string model_file_name = argv[1];
    std::string output_file_name = argv[2];
    std::string acc_structure = argv[3];

    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int max_depth = 50;

    // World - Read our .obj file
    std::vector<std::unique_ptr<const Mesh>> meshes;
    makeScene(meshes, model_file_name);

    std::unique_ptr<AccelerationStructure> accel;
    if (acc_structure == "none")
    {
        accel = std::unique_ptr<AccelerationStructure>(new AccelerationStructure(meshes));
    }
    else if (acc_structure == "grid")
    {
        accel = std::unique_ptr<AccelerationStructure>(new Grid(meshes));
    }
    else if (acc_structure == "bvh")
    {
        // accel = std::unique_ptr<AccelerationStructure>(new BVH(meshes));
        std::cerr << "Not implemented yet" << std::endl;
        return 2;
    }
    else
    {
        std::cerr << "Not implemented yet" << std::endl;
        return 2;
    }

    // Camera
    point3 lookfrom(13, 2, 3);
    point3 lookat(0, 0, 0);
    Vec3f vup(0, 1, 0);
    auto dist_to_focus = 10.0;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, dist_to_focus);

    // Render - store to PPM file
    std::ofstream ofs;
    ofs.open(output_file_name);
    ofs << "P3\n"
        << image_width << ' ' << image_height << "\n255\n";

    const clock_t begin_time = clock();
    for (int j = image_height - 1; j >= 0; --j)
    {
        std::cout << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i)
        {
            auto u = double(i) / (image_width - 1);
            auto v = double(j) / (image_height - 1);

            ray r = cam.get_ray(u, v);

            color pixel_color = ray_color(r, accel, cam, max_depth);

            write_color(ofs, pixel_color);
        }
    }
    time_diff = float(clock() - begin_time) / CLOCKS_PER_SEC;
    ofs.close();
    cout.precision(2);
    std::cout << "\nDone.\nCzas: " << time_diff << "[s]\n";
}
