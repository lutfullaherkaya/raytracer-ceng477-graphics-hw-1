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

#define MAX_DEPTH 19

using namespace parser;


struct BVHNode {
    BVHNode() = default;

    int depth = 0;
    std::list<Triangle> faces; // is empty if not leaf.
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
            node->faces = faces; // is leaf node
            node->spheres = spheres;
        } else {
            std::list<Triangle> leftHalf = {};
            std::list<Triangle> rightHalf = {};
            std::list<Sphere> leftSprs = {};
            std::list<Sphere> rightSprs = {};

            int axis = depth % 3;
            Box currentBox = node->box;
            auto start = currentBox.min[axis], end = currentBox.max[axis];
            auto midPoint = (start + end) / 2;
            int maxTries = 10; // doing this to eliminate empty boxes
            while (maxTries-- && ((leftHalf.empty() && leftSprs.empty()) || (rightHalf.empty() && rightSprs.empty()))) {
                leftHalf.clear();
                rightHalf.clear();
                leftSprs.clear();
                rightSprs.clear();

                for (auto &face: faces) {// splitting by location instad of sorting the faces and splitting from median is way faster
                    auto vertexPoint = scene.getCenter(face.indices)[axis];
                    if (vertexPoint < midPoint) {
                        leftHalf.push_back(face);
                    } else {
                        rightHalf.push_back(face);
                    }
                }
                for (auto &sphere: spheres) {// splitting by location instad of sorting the faces and splitting from median is way faster
                    auto vertexPoint = scene.vertex_data[sphere.center_vertex_id - 1][axis];
                    if (vertexPoint < midPoint) {
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

            if ((leftHalf.empty() && leftSprs.empty()) || (rightHalf.empty() && rightSprs.empty())) {
                node->faces = faces; // is leaf node now because we couldn't split it by space
            } else {
                node->right = build(rightHalf, rightSprs, depth + 1, scene);
                node->left = build(leftHalf, leftSprs, depth + 1, scene);
            }

        }

        return node;
    }



};


struct BVHTree {
    BVHNode *root;
    Scene &scene;

    explicit BVHTree(Scene &scene) : root(nullptr), scene(scene) {}

    void build(std::list<Triangle> triangles) {
        root = BVHNode::build(std::move(triangles), {scene.spheres.begin(), scene.spheres.end()}, 0, scene);
    }



};


#endif //RAYTRACER_CENG477_GRAPHICS_HW_1_BVH_H
