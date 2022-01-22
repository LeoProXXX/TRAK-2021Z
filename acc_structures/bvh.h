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

    struct Tree
    {
        Tree(const Extents& sceneExtents, uint32_t& max_tree_depth, uint32_t& children_num);
        ~Tree();
        void insert(const Extents* extents);
        void build();

        struct TreeNode
        {
            TreeNode(uint32_t* children_num);
            TreeNode* child[(*children_num)] = { nullptr };;
            std::vector<const Extents*> nodeExtentsList;
            Extents nodeExtents;
            bool isLeaf = true;
            uint32_t* children_n_adr;
        };

        struct QueueElement
        {
            const TreeNode* node; 
            float t; 
            QueueElement(const TreeNode* n, float tn) : node(n), t(tn) {}
            friend bool operator < (const QueueElement& a, const QueueElement& b) { return a.t > b.t; }
        };
        
        TreeNode* root = nullptr;
        BBox<> bbox;
        uint32_t tree_depth;
        uint32_t children_no;

    private:

        void deleteTreeNode(TreeNode*& node);
        void insert(TreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth);
        void build(TreeNode*& node, const BBox<>& bbox);
    };

    std::vector<Extents> extentsList;
    Tree* tree = nullptr;

public:
    BVH(std::vector<std::unique_ptr<const Mesh>>& m, uint32_t max_tree_depth, uint32_t children_num);
    bool intersect(const Vec3f&, const Vec3f&, const uint32_t&, float&, hit_record&) const;
    ~BVH();
};

#endif