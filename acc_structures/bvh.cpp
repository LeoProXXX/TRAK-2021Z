#include "bvh.h"

const Vec3f BVH::planeSetNormals[BVH::kNumPlaneSetNormals] = {
    Vec3f(1, 0, 0),
    Vec3f(0, 1, 0),
    Vec3f(0, 0, 1),
    Vec3f(sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3f(-sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3f(-sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f),
    Vec3f(sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f)
};


BVH::BVH(std::vector<std::unique_ptr<const Mesh>>& m, uint32_t max_tree_depth, uint32_t children_num) : AccelerationStructure(m)
{
    Extents sceneExtents;
    extentsList.reserve(meshes.size());
    for (uint32_t i = 0; i < meshes.size(); ++i) {
        for (uint8_t j = 0; j < kNumPlaneSetNormals; ++j) {
            for (const auto vtx : meshes[i]->vertexPool) {
                float d = dot(planeSetNormals[j], vtx);

                if (d < extentsList[i].d[j][0]) extentsList[i].d[j][0] = d;
                if (d > extentsList[i].d[j][1]) extentsList[i].d[j][1] = d;
            }
        }
        sceneExtents.extendBy(extentsList[i]);
        extentsList[i].mesh = meshes[i].get();
    }

    tree = new Tree(sceneExtents, max_tree_depth, children_num);

    for (uint32_t i = 0; i < meshes.size(); ++i) {
        tree->insert(&extentsList[i]);
    }

    tree->build();
}

BVH::Extents::Extents()
{
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i)
        d[i][0] = kInfinity, d[i][1] = -kInfinity;
}

void BVH::Extents::extendBy(const Extents& e)
{

    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        if (e.d[i][0] < d[i][0]) d[i][0] = e.d[i][0];
        if (e.d[i][1] > d[i][1]) d[i][1] = e.d[i][1];
    }
}

Vec3f BVH::Extents::centroid() const
{
    return Vec3f(
        d[0][0] + d[0][1] * 0.5,
        d[1][0] + d[1][1] * 0.5,
        d[2][0] + d[2][1] * 0.5);
}

BVH::Tree::Tree(const Extents& sceneExtents, uint32_t& max_tree_depth, uint32_t& children_num)
{
    float xDiff = sceneExtents.d[0][1] - sceneExtents.d[0][0];
    float yDiff = sceneExtents.d[1][1] - sceneExtents.d[1][0];
    float zDiff = sceneExtents.d[2][1] - sceneExtents.d[2][0];
    float maxDiff = std::max(xDiff, std::max(yDiff, zDiff));
    Vec3f minPlusMax(
        sceneExtents.d[0][0] + sceneExtents.d[0][1],
        sceneExtents.d[1][0] + sceneExtents.d[1][1],
        sceneExtents.d[2][0] + sceneExtents.d[2][1]);
    bbox[0] = (minPlusMax - maxDiff) * 0.5;
    bbox[1] = (minPlusMax + maxDiff) * 0.5;
    this->tree_depth = max_tree_depth;
    this->children_no = children_num;
    root = new TreeNode(&this->children_no);
}

BVH::Tree::~Tree() { deleteTreeNode(root); }

void BVH::Tree::insert(const Extents* extents) { insert(root, extents, bbox, 0); }

void BVH::Tree::build() { build(root, bbox); };

BVH::Tree::TreeNode::TreeNode(uint32_t* children_num) 
{ 
this->children_n_adr = children_num;
std::cout << "childern numb: " << *children_num << std::endl;
std::cout << "childern numb: " << *(this->children_n_adr) << std::endl;
}

void BVH::Tree::deleteTreeNode(TreeNode*& node)
{
    for (uint8_t i = 0; i < 8; i++) {
        if (node->child[i] != nullptr) {
            deleteTreeNode(node->child[i]);
        }
    }
    delete node;
}

void BVH::Tree::insert(TreeNode*& node, const Extents* extents, const BBox<>& bbox, uint32_t depth)
{
    if (node->isLeaf) {
        if (node->nodeExtentsList.size() == 0 || depth == this->tree_depth) {
            node->nodeExtentsList.push_back(extents);
        }
        else {
            node->isLeaf = false;

            while (node->nodeExtentsList.size()) {
                insert(node, node->nodeExtentsList.back(), bbox, depth);
                node->nodeExtentsList.pop_back();
            }

            insert(node, extents, bbox, depth);
        }
    }
    else {
        Vec3f extentsCentroid = extents->centroid();
        Vec3f nodeCentroid = (bbox[0] + bbox[1]) * 0.5;
        BBox<> childBBox;
        uint8_t childIndex = 0;
        // x-axis
        if (extentsCentroid.x > nodeCentroid.x) {
            childIndex = 4;
            childBBox[0].x = nodeCentroid.x;
            childBBox[1].x = bbox[1].x;
        }
        else {
            childBBox[0].x = bbox[0].x;
            childBBox[1].x = nodeCentroid.x;
        }
        // y-axis
        if (extentsCentroid.y > nodeCentroid.y) {
            childIndex += 2;
            childBBox[0].y = nodeCentroid.y;
            childBBox[1].y = bbox[1].y;
        }
        else {
            childBBox[0].y = bbox[0].y;
            childBBox[1].y = nodeCentroid.y;
        }
        // z-axis
        if (extentsCentroid.z > nodeCentroid.z) {
            childIndex += 1;
            childBBox[0].z = nodeCentroid.z;
            childBBox[1].z = bbox[1].z;
        }
        else {
            childBBox[0].z = bbox[0].z;
            childBBox[1].z = nodeCentroid.z;
        }

        if (node->child[childIndex] == nullptr)
            node->child[childIndex] = new TreeNode(&this->children_no);
        insert(node->child[childIndex], extents, childBBox, depth + 1);
    }
}

void BVH::Tree::build(TreeNode*& node, const BBox<>& bbox)
{
    if (node->isLeaf) {
        for (const auto& e : node->nodeExtentsList) {
            node->nodeExtents.extendBy(*e);
        }
    }
    else {
        for (uint8_t i = 0; i < 8; ++i) {
            if (node->child[i]) {
                BBox<> childBBox;
                Vec3f centroid = bbox.centroid();
                // x-axis
                childBBox[0].x = (i & 4) ? centroid.x : bbox[0].x;
                childBBox[1].x = (i & 4) ? bbox[1].x : centroid.x;
                // y-axis
                childBBox[0].y = (i & 2) ? centroid.y : bbox[0].y;
                childBBox[1].y = (i & 2) ? bbox[1].y : centroid.y;
                // z-axis
                childBBox[0].z = (i & 1) ? centroid.z : bbox[0].z;
                childBBox[1].z = (i & 1) ? bbox[1].z : centroid.z;

                // Inspect child
                build(node->child[i], childBBox);

                node->nodeExtents.extendBy(node->child[i]->nodeExtents);
            }
        }
    }
}

bool BVH::Extents::intersect(
    const float* precomputedNumerator,
    const float* precomputedDenominator,
    float& tNear,
    float& tFar,
    uint8_t& planeIndex) const
{
    numRayBoundingVolumeTests++;
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        float tNearExtents = (d[i][0] - precomputedNumerator[i]) / precomputedDenominator[i];
        float tFarExtents = (d[i][1] - precomputedNumerator[i]) / precomputedDenominator[i];
        if (precomputedDenominator[i] < 0) std::swap(tNearExtents, tFarExtents);
        if (tNearExtents > tNear) tNear = tNearExtents, planeIndex = i;
        if (tFarExtents < tFar) tFar = tFarExtents;
        if (tNear > tFar) return false;
    }

    return true;
}

bool BVH::intersect(const Vec3f& orig, const Vec3f& dir, const uint32_t& rayId, float& tHit, hit_record& rec) const
{
    tHit = kInfinity;
    const Mesh* intersectedMesh = nullptr;
    float precomputedNumerator[BVH::kNumPlaneSetNormals];
    float precomputedDenominator[BVH::kNumPlaneSetNormals];
    for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
        precomputedNumerator[i] = dot(planeSetNormals[i], orig);
        precomputedDenominator[i] = dot(planeSetNormals[i], dir);
    }

    uint8_t planeIndex;
    float tNear = 0, tFar = kInfinity; 
    if (!tree->root->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNear, tFar, planeIndex) || tFar < 0)
        return false;
    tHit = tFar;
    std::priority_queue<BVH::Tree::QueueElement> queue;
    queue.push(BVH::Tree::QueueElement(tree->root, 0));
    while (!queue.empty() && queue.top().t < tHit) {
        const Tree::TreeNode* node = queue.top().node;
        queue.pop();
        if (node->isLeaf) {
            for (const auto& e : node->nodeExtentsList) {
                float t = kInfinity;
                if (e->mesh->intersect(orig, dir, t, rec) && t < tHit) {
                    tHit = t;
                    intersectedMesh = e->mesh;
                }
            }
        }
        else {
            for (uint8_t i = 0; i < 8; ++i) {
                if (node->child[i] != nullptr) {
                    float tNearChild = 0, tFarChild = tFar;
                    if (node->child[i]->nodeExtents.intersect(precomputedNumerator, precomputedDenominator, tNearChild, tFarChild, planeIndex)) {
                        float t = (tNearChild < 0 && tFarChild >= 0) ? tFarChild : tNearChild;
                        queue.push(BVH::Tree::QueueElement(node->child[i], t));
                    }
                }
            }
        }
    }

    return (intersectedMesh != nullptr);
}

BVH::~BVH() { delete tree; }