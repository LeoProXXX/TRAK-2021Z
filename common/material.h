#ifndef MATERIAL_H
#define MATERIAL_H

#include "helpers.h"

struct hit_record;

class material
{

public:
    material(Vec3f ambient_param, Vec3f diffuse_param, Vec3f specular_param, float shininess_param)
        : ambient(ambient_param), diffuse(diffuse_param), specular(specular_param), shininess(shininess_param){};

public:
    Vec3f ambient;
    Vec3f diffuse;
    Vec3f specular;
    float shininess;
};

#endif
