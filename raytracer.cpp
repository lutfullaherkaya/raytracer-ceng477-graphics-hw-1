#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <chrono>
#include <list>
#include <stack>

/* optimization ideas:
 * - dont check each mesh one by one, instead use another tree structure to store the meshes
 * - use bounding boxes instead of bounding spheres in tree structure
 * */


using namespace parser;


float det(float m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

Vec3f getCenter(Face &face, Scene &scene) {
    Vec3f center = scene.vertex_data[face.v0_id-1] + scene.vertex_data[face.v1_id-1] + scene.vertex_data[face.v2_id-1];
    return center / 3;
}

struct IntersectionPoint {
    float t1, t2;
    //Vec3f normal;
    bool exists;
};

class Ray {
    Vec3f getPoint(float t) { // t >= 0
        return origin + direction * t;
    }

    Vec3f origin;
    Vec3f direction;

public:
    Ray(Vec3f origin, Vec3f direction) : origin(origin), direction(direction) {

    }


    IntersectionPoint intersects(Scene &scene, MySphere mySphere) {
        auto c = mySphere.c;
        auto r = mySphere.r;
        auto d = direction;
        auto o = origin;

        auto B = 2 * (d * (o - c));
        auto A = d * d;
        auto C = (o - c) * (o - c) - r * r;
        auto discriminant = B * B - 4 * A * C;

        if (discriminant >= 0) {
            float t1 = (-d * (o - c) - sqrt(discriminant)) / (2 * (d * d));
            float t2 = (-d * (o - c) + sqrt(discriminant)) / (2 * (d * d));
            return {t1, t2, true};
        }
        return {0, 0, false};

    }

    // source: book
    IntersectionPoint intersects(Scene &scene, Box box) {
        float tmin = (box.min.x - origin.x) / direction.x;
        float tmax = (box.max.x - origin.x) / direction.x;

        if (tmin > tmax) std::swap(tmin, tmax);

        float tymin = (box.min.y - origin.y) / direction.y;
        float tymax = (box.max.y - origin.y) / direction.y;

        if (tymin > tymax) std::swap(tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax))
            return {0, 0, false};

        if (tymin > tmin)
            tmin = tymin;

        if (tymax < tmax)
            tmax = tymax;

        float tzmin = (box.min.z - origin.z) / direction.z;
        float tzmax = (box.max.z - origin.z) / direction.z;

        if (tzmin > tzmax) std::swap(tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax))
            return {0, 0, false};

        if (tzmin > tmin)
            tmin = tzmin;

        if (tzmax < tmax)
            tmax = tzmax;

        return {tmin, tmax, true};

    }

    IntersectionPoint intersects(Scene &scene, Face &face, int material_id) {
        IntersectionPoint point = {0, 0, false};
        auto a = scene.vertex_data[face.v0_id - 1];
        auto b = scene.vertex_data[face.v1_id - 1];
        auto c = scene.vertex_data[face.v2_id - 1];

        float A[3][3] = {
                {a.x - b.x, a.x - c.x, direction.x},
                {a.y - b.y, a.y - c.y, direction.y},
                {a.z - b.z, a.z - c.z, direction.z}
        };
        float detA = det(A);

        float betaMatrix[3][3] = {
                {a.x - origin.x, a.x - c.x, direction.x},
                {a.y - origin.y, a.y - c.y, direction.y},
                {a.z - origin.z, a.z - c.z, direction.z}
        };
        float beta = det(betaMatrix) / detA;

        float gammaMatrix[3][3] = {
                {a.x - b.x, a.x - origin.x, direction.x},
                {a.y - b.y, a.y - origin.y, direction.y},
                {a.z - b.z, a.z - origin.z, direction.z}
        };
        float gamma = det(gammaMatrix) / detA;

        float tMatrix[3][3] = {
                {a.x - b.x, a.x - c.x, a.x - origin.x},
                {a.y - b.y, a.y - c.y, a.y - origin.y},
                {a.z - b.z, a.z - c.z, a.z - origin.z}
        };
        float t = det(tMatrix) / detA;

        float alpha = 1 - beta - gamma;
        float tMin = 0;

        if (alpha >= 0 &&
            beta >= 0 &&
            gamma >= 0 &&
            t >= tMin) {
            point.exists = true;
            point.t1 = t;
        }
        return point;

    }
};

#define MAX_DEPTH 19


struct BVHNode {
    BVHNode() = default;

    int depth = 0;
    std::list<Face> faces; // is empty if not leaf.
    Box box{};
    BVHNode *left = nullptr, *right = nullptr;

    static BVHNode *build(std::list<Face> faces, int depth, Scene &scene) {
        if (faces.empty()) {
            return nullptr;
        }

        auto *node = new BVHNode();
        node->depth = depth;
        node->box = scene.getBoundingBox(faces);

        if (faces.size() <= 1 || depth >= MAX_DEPTH) {
            node->faces = faces; // is leaf node
        } else {
            std::list<Face> leftHalf = {};
            std::list<Face> rightHalf = {};

            int axis = depth % 3;


            Box currentBox = node->box;
            auto start = currentBox.min[axis], end = currentBox.max[axis];
            auto midPoint = (start + end) / 2;
            int maxTries = 10; // doing this to eliminate empty boxes
            while (maxTries-- && (leftHalf.empty() || rightHalf.empty())) {
                leftHalf.clear();
                rightHalf.clear();

                for (auto &face: faces) {// splitting by location instad of sorting the faces and splitting from median is way faster
                    auto vertexPoint = getCenter(face, scene)[axis];
                    if (vertexPoint < midPoint) {
                        leftHalf.push_back(face);
                    } else {
                        rightHalf.push_back(face);
                    }
                }

                if (leftHalf.empty()) {
                    start = midPoint;
                    midPoint = (start + end) / 2;
                }
                if (rightHalf.empty()) {
                    end = midPoint;
                    midPoint = (start + end) / 2;
                }
            }

            if (leftHalf.empty() || rightHalf.empty()) {
                node->faces = faces; // is leaf node now because we couldn't split it by space
            } else {
                node->right = build(rightHalf, depth + 1, scene);
                node->left = build(leftHalf, depth + 1, scene);
            }

        }

        return node;
    }


};


struct BVHTree {
    BVHNode *root;
    Mesh &mesh;
    Scene &scene;

    BVHTree(Mesh &mesh, Scene &scene) : mesh(mesh), scene(scene) {
        std::list<Face> faceList(mesh.faces.begin(), mesh.faces.end());
        root = BVHNode::build(faceList, 0, scene);
    }

};

class RayTracer {
public:
    RayTracer(parser::Scene &scene) {
        this->scene = scene;
    }

    std::vector<BVHTree> meshTrees;

    void rayTrace() {
        for (auto camera: scene.cameras) {
            auto *image = new unsigned char[camera.image_width * camera.image_height * 3];
            int imagePtr = 0;
            for (int j = 0; j < camera.image_height; j++) {
                std::cout << 100 * j / (double) camera.image_height << '%' << std::endl;
                for (int i = 0; i < camera.image_width; i++) {
                    Ray eyeRay = generateEyeRay(camera, i, j);

                    auto raytracedColor = scene.background_color;

                    for (auto &sphere: scene.spheres) {
                        auto intersectionPoint = eyeRay.intersects(scene,
                                                                   {scene.vertex_data[sphere.center_vertex_id - 1],
                                                                    sphere.radius});
                        if (intersectionPoint.exists) {
                            raytracedColor.x = 255;
                            raytracedColor.y = 0;
                            raytracedColor.z = 0;
                        }
                    }
                    for (auto &triangle: scene.triangles) {
                        auto intersectionPoint = eyeRay.intersects(scene, triangle.indices, triangle.material_id);
                        if (intersectionPoint.exists) {
                            raytracedColor.x = 255;
                            raytracedColor.y = 255;
                            raytracedColor.z = 0;
                        }
                    }
                    for (int meshIndex = 0; meshIndex < scene.meshes.size(); meshIndex++) {

                        std::stack<BVHNode *> stack;
                        stack.push(meshTrees[meshIndex].root);
                        while (!stack.empty()) {
                            auto node = stack.top();
                            stack.pop();
                            if (eyeRay.intersects(scene, node->box).exists) {
                                if (node->left) {
                                    stack.push(node->left);
                                }
                                if (node->right) {
                                    stack.push(node->right);
                                }
                                if (!node->left && !node->right) {
                                    for (auto face: node->faces) {
                                        auto &mesh = scene.meshes[meshIndex];
                                        auto intersectionPoint = eyeRay.intersects(scene, face, mesh.material_id);
                                        if (intersectionPoint.exists) {
                                            raytracedColor.x = 0;
                                            raytracedColor.y = 0;
                                            raytracedColor.z = 255;
                                        }
                                    }
                                }
                            }
                        }


                    }

                    image[imagePtr++] = raytracedColor.x;
                    image[imagePtr++] = raytracedColor.y;
                    image[imagePtr++] = raytracedColor.z;
                }
            }
            write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);
            std::cout << camera.image_name << std::endl;

        }

    }

    static Ray generateEyeRay(Camera &camera, int i, int j) {
        auto e = camera.position;
        auto w = -camera.gaze;
        auto distance = camera.near_distance;

        auto l = camera.near_plane.x;
        auto r = camera.near_plane.y;
        auto b = camera.near_plane.z;
        auto t = camera.near_plane.w;

        auto v = camera.up;
        auto u = v.crossProduct(w);
        auto nx = camera.image_width;
        auto ny = camera.image_height;

        auto m = e + -w * distance;
        auto q = m + u * l + v * t;

        float su = (i + 0.5) * (r - l) / nx;
        float sv = (j + 0.5) * (t - b) / ny;

        auto s = q + u * su - v * sv;

        return {e, s - e};
    }

private:
    parser::Scene scene;

};

int main(int argc, char *argv[]) {
    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    RayTracer raytracer(scene);


    auto begin1 = std::chrono::high_resolution_clock::now();
    for (auto &mesh: scene.meshes) {
        raytracer.meshTrees.emplace_back(mesh, scene);
    }
    auto end1 = std::chrono::high_resolution_clock::now();

    auto begin2 = std::chrono::high_resolution_clock::now();
    raytracer.rayTrace();
    auto end2 = std::chrono::high_resolution_clock::now();

    auto elapsed1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1);
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2);
    printf("\nPlanted trees in %.3f seconds.\n", elapsed1.count() * 1e-9);
    printf("Rendered in %.3f seconds.\n", elapsed2.count() * 1e-9);
}
