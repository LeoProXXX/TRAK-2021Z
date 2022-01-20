#ifndef BVH_H
#define BVH_H

#include "../common/helpers.h"
#include "../common/mesh.h"
#include "acceleration_structure.h"

class BVH : public AccelerationStructure
{
	static const uint8_t kNumPlaneSetNormals = 7;
	static const Vec3f planeSetNormals[kNumPlaneSetNormals];
    struct Extents
    {
        Extents();
        void extendBy(const Extents& e);
        Vec3f centroid() const;
        bool intersect(const float*, const float*, float&, float&, uint8_t&) const;
        float d[kNumPlaneSetNormals][2];
        const Mesh* mesh;
    };

    struct Octree
    {
        Octree(const Extents& sceneExtents);
        ~Octree();
        void insert(const Extents* extents);
        void build();

        struct OctreeNode
        {
            OctreeNode* child[8] = { nullptr };;
            std::vector<const Extents*> nodeExtentsList;
            Extents nodeExtents;
            bool isLeaf = true;
        };

        struct QueueElement
        {
            const OctreeNode* node; 
            float t; 
            QueueElement(const OctreeNode* n, float tn) : node(n), t(tn) {}

            friend bool operator < (const QueueElement& a, const QueueElement& b) { return a.t > b.t; }
        };
        
        OctreeNode* root = nullptr;
        BBox<> bbox;

    private:

        void deleteOctreeNode(OctreeNode*& node);
        void insert(OctreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth);
        void build(OctreeNode*& node, const BBox<>& bbox);
    };

    std::vector<Extents> extentsList;
    Octree* octree = nullptr;

public:
	BVH(std::vector<std::unique_ptr<const Mesh>>& m);
    bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&, hit_record&) const;
    ~BVH();
};

#endif