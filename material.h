#ifndef MATERIAL_H
#define MATERIAL_H

#include "common/helpers.h"

struct hit_record;

class material
{

public:
    material(vec3 ambient_param, vec3 diffuse_param, vec3 specular_param, float shininess_param)
        : ambient(ambient_param), diffuse(diffuse_param), specular(specular_param), shininess(shininess_param){};

public:
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float shininess;
};

#endif
