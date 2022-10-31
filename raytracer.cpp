#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <chrono>
#include <list>
#include <stack>


using namespace parser;


float det(float m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
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

#define MAX_DEPTH 8

struct SphereTreeNode {
    SphereTreeNode() = default;

    int depth = 0;
    std::list<Face> faces; // is empty if not leaf.
    MySphere sphere{};
    SphereTreeNode *left = nullptr, *right = nullptr;

    static SphereTreeNode *build(std::list<Face> faces, int depth, Scene &scene) {
        if (faces.empty()) {
            return nullptr;
        }

        auto *node = new SphereTreeNode();
        node->depth = depth;
        node->sphere = scene.getBoundingSphere(faces);

        if (faces.size() <= 1 || depth >= MAX_DEPTH) {
            node->faces = faces; // is leaf node
        } else {

            int axis = depth % 3;
            // sort faces
            faces.sort([scene, axis](Face &a, Face &b) {
                auto vertexA = scene.vertex_data[a.v0_id - 1];
                auto vertexB = scene.vertex_data[b.v0_id - 1];
                return vertexA[axis] < vertexB[axis];
            });

            std::list<Face> &leftHalf = faces;
            std::list<Face> rightHalf = {};
            auto halfItr = std::next(leftHalf.begin(), leftHalf.size() / 2);
            leftHalf.splice(rightHalf.begin(), leftHalf, halfItr, leftHalf.end());

            node->right = build(rightHalf, depth + 1, scene);
            node->left = build(leftHalf, depth + 1, scene);

        }

        return node;
    }


};


struct SphereTree {
    SphereTreeNode *root;
    Mesh &mesh;
    Scene &scene;

    SphereTree(Mesh &mesh, Scene &scene) : mesh(mesh), scene(scene) {
        std::list<Face> faceList(mesh.faces.begin(), mesh.faces.end());
        root = SphereTreeNode::build(faceList, 0, scene);
    }

};

class RayTracer {
public:
    RayTracer(parser::Scene &scene) {
        this->scene = scene;
    }

    std::vector<SphereTree> meshSphereTrees;

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

                        std::stack<SphereTreeNode *> stack;
                        stack.push(meshSphereTrees[meshIndex].root);
                        while (!stack.empty()) {
                            auto node = stack.top();
                            stack.pop();
                            if (eyeRay.intersects(scene, node->sphere).exists) {
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
        raytracer.meshSphereTrees.emplace_back(mesh, scene);
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
