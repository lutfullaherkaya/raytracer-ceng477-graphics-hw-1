//
// Created by lutfu on 1.11.2022.
//

#ifndef RAYTRACER_CENG477_GRAPHICS_HW_1_BVH_H
#define RAYTRACER_CENG477_GRAPHICS_HW_1_BVH_H

#include <cmath>
#include <chrono>
#include <list>
#include <stack>
#include <utility>
#include "parser.h"
#include <thread>
#include <iostream>

#define MAX_DEPTH 19

using namespace parser;

/**
 * How I implemented this? - Lütfullah Erkaya
 * At first I tried implementing a KD tree but it didn't work. I was dividing the space with points and not triangles.
 * Then I thought dividing the triangles, not points was simpler. Only thing I had to do was dividing the triangles list
 * into two and then recursively do the same thing for the two lists. I also had to calculate the bounding box of the
 * triangles list. Then the tree was formed. I thought it was still a KD tree but it was BVH.
 *
 * For dividing in half, at first i sorted the list and divided by median but it was very slow, it was taking 4 seconds
 * to build the tree for horse_and_mug. Then I read in course book that it was faster to not sort but just partition using
 * the middle point of the parent box.
 * Time it took to plant went from 4 seconds to 0,08 seconds.
 *
 * Then I did the divide by widest axis and first box intersection traversal, the rendering became much faster.
 *
 */
struct BVHNode {
    BVHNode() = default;

    int depth = 0;
    int axis = 0;
    std::list<Triangle> triangles; // is empty if not leaf.
    std::list<Sphere> spheres; // is empty if not leaf.
    Box box{};
    BVHNode *left = nullptr, *right = nullptr;

    static BVHNode *build(std::list<Triangle> faces, std::list<Sphere> spheres, int depth, Scene &scene) {
        if (faces.empty() && spheres.empty()) {
            return nullptr;
        }
        auto *node = new BVHNode();
        node->depth = depth;
        node->box = scene.getBoundingBox(faces);
        scene.extendBoundingBox(node->box, spheres);

        if (faces.size() + spheres.size() <= 1 || depth >= MAX_DEPTH) {
            node->triangles = faces; // is leaf node
            node->spheres = spheres;
        } else {
            node->axis = node->box.getWidestAxis();

            std::list<Triangle> leftHalf = {};
            std::list<Triangle> rightHalf = {};
            std::list<Sphere> leftSprs = {};
            std::list<Sphere> rightSprs = {};

            if (partition(node, faces, spheres, leftHalf, leftSprs, rightHalf, rightSprs, scene)) {
                node->right = build(rightHalf, rightSprs, depth + 1, scene);
                node->left = build(leftHalf, leftSprs, depth + 1, scene);
            } else {
                node->triangles = faces; // is leaf node now because we couldn't split it by space
                node->spheres = spheres;
            }

        }

        return node;
    }

    bool isLeaf() {
        return left == nullptr && right == nullptr;
    }

    static bool
    partition(BVHNode *node, std::list<Triangle> &triangles, std::list<Sphere> &spheres, std::list<Triangle> &leftHalf,
              std::list<Sphere> &leftSprs, std::list<Triangle> &rightHalf, std::list<Sphere> &rightSprs, Scene &scene) {
        Box currentBox = node->box;
        auto start = currentBox.min[node->axis], end = currentBox.max[node->axis];
        auto midPoint = (start + end) / 2;
        int maxTries = 19; // doing this to eliminate empty boxes. this is very important for horse_and_mug. withot this, it does not render.
        while (maxTries-- && ((leftHalf.empty() && leftSprs.empty()) || (rightHalf.empty() && rightSprs.empty()))) {
            leftHalf.clear();
            rightHalf.clear();
            leftSprs.clear();
            rightSprs.clear();

            for (auto &face: triangles) {// splitting by location instad of sorting the triangles and splitting from median is way faster
                if (face.center[node->axis] < midPoint) {
                    leftHalf.push_back(face);
                } else {
                    rightHalf.push_back(face);
                }
            }
            for (auto &sphere: spheres) {// splitting by location instad of sorting the triangles and splitting from median is way faster
                if (scene.vertex_data[sphere.center_vertex_id - 1][node->axis] < midPoint) {
                    leftSprs.push_back(sphere);
                } else {
                    rightSprs.push_back(sphere);
                }
            }

            if (leftHalf.empty() && leftSprs.empty()) {
                start = midPoint;
                midPoint = (start + end) / 2;
            }
            if (rightHalf.empty() && rightSprs.empty()) {
                end = midPoint;
                midPoint = (start + end) / 2;
            }
        }
        return !((leftHalf.empty() && leftSprs.empty()) || (rightHalf.empty() && rightSprs.empty()));
    }
};


struct BVHTree {
    BVHNode *root;
    Scene &scene;

    explicit BVHTree(Scene &scene) : root(nullptr), scene(scene) {}

    void build(std::list<Triangle> triangles) {
        root = BVHNode::build(std::move(triangles),
                              {scene.spheres.begin(), scene.spheres.end()}, 0, scene);
    }


};


#endif //RAYTRACER_CENG477_GRAPHICS_HW_1_BVH_H
