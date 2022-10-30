#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <chrono>
#include <algorithm>
#include <list>
#include <utility>

using namespace parser;


float det(float m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

float min(float a, float b) {
    return a < b ? a : b;
}

float max(float a, float b) {
    return a > b ? a : b;
}

void swap(float &a, float &b) {
    float temp = a;
    a = b;
    b = temp;
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

    bool intersects(Box box) const {
        float txmin = (box.min.x - origin.x) / direction.x;
        float txmax = (box.max.x - origin.x) / direction.x;
        if (txmin > txmax) swap(txmin, txmax);

        float tymin = (box.min.y - origin.y) / direction.y;
        float tymax = (box.max.y - origin.y) / direction.y;
        if (tymin > tymax) swap(tymin, tymax);

        if ((txmin > tymax) || (tymin > txmax))
            return false;

        if (tymin > txmin)
            txmin = tymin;

        if (tymax < txmax)
            txmax = tymax;

        float tzmin = (box.min.z - origin.z) / direction.z;
        float tzmax = (box.max.z - origin.z) / direction.z;
        if (tzmin > tzmax) swap(tzmin, tzmax);

        if ((txmin > tzmax) || (tzmin > txmax))
            return false;

        return true;
    }
};


class KdThree3dNode {
public:
    KdThree3dNode *left, *right;
    Box box;
    int depth;
    static int maxDepth;
    std::list<int> tVertexIndices;

    KdThree3dNode() : left(nullptr), right(nullptr), depth(0) {

    }

    KdThree3dNode(Box box, int depth, int startIndex, int endIndex, std::list<int> tVertexIndices)
            : left(nullptr), right(nullptr), box(box), depth(depth), tVertexIndices(std::move(tVertexIndices)) {

    }

    static KdThree3dNode *build(int depth, std::list<int> tVertexIndices, Scene &scene,
                                std::vector<TriangleVertex> &tVertices,
                                Box box,
                                bool isLeftHalf,
                                float midPoint) {
        auto *node = new KdThree3dNode();
        node->depth = depth;
        node->tVertexIndices = tVertexIndices;
        node->box = box;

        int axis = depth % 3;

        if (node->tVertexIndices.size() <= 1) {
            return node;
        }

        if (depth + 1 < maxDepth) {
            int newAxis = (depth + 1) % 3;
            int median = (node->tVertexIndices.size() - 1) / 2;

            node->tVertexIndices.sort([&tVertices, newAxis, scene](int a, int b) {
                auto aVertex = scene.vertex_data[tVertices[a].v_id - 1];
                auto bVertex = scene.vertex_data[tVertices[b].v_id - 1];
                return aVertex[newAxis] < bVertex[newAxis];
            });
            std::list<int> leftHalf = tVertexIndices;
            std::list<int> rightHalf;
            leftHalf.splice(rightHalf.begin(), leftHalf, std::next(leftHalf.begin(), median + 1), leftHalf.end());


            // check if box is 3d
            Box leftBox = box;
            Box rightBox = box;
            auto mediumPoint = scene.vertex_data[tVertices[leftHalf.back()].v_id - 1][newAxis];
            leftBox.max[newAxis] = mediumPoint;
            rightBox.min[newAxis] = mediumPoint;

            if (leftBox.dimensions() > Vec3f(0.01, 0.01, 0.01) && rightBox.dimensions() > Vec3f(0.01, 0.01, 0.01) &&
                leftBox.is3d() && rightBox.is3d()) {
                node->left = build(depth + 1, leftHalf, scene, tVertices,
                                   leftBox, true, mediumPoint);
                node->right = build(depth + 1, rightHalf, scene, tVertices,
                                    rightBox, false, mediumPoint);
            }


        }


        return node;


    }

    std::vector<KdThree3dNode *> intersectingNodes(Ray &ray) {
        if (!ray.intersects(box)) {
            return {};
        }

        if (isLeaf()) {
            return {this};
        }

        std::vector<KdThree3dNode *> leftNodes;
        std::vector<KdThree3dNode *> rightNodes;
        if (left) {
            leftNodes = left->intersectingNodes(ray);
        }

        if (right) {
            rightNodes = right->intersectingNodes(ray);
        }

        leftNodes.insert(leftNodes.end(), rightNodes.begin(), rightNodes.end());
        return leftNodes;
    }

    ~KdThree3dNode() {
        delete left;
        delete right;
    }

    bool isLeaf() {
        return left == nullptr && right == nullptr;
    }
};


class OctaTreeNode {
public:
    OctaTreeNode *nextNode;
    Box box;
    int depth;
    std::vector<int> v_ids;

    static OctaTreeNode *build(std::vector<int> v_ids_new, Scene &scene, Box box, int depth) {
        if (v_ids_new.size() > 0) {
            auto *node = new OctaTreeNode();
            node->box = box;
            node->depth = depth;
            node->nextNode = nullptr;
            node->v_ids = v_ids_new;

            if (depth < 3) {
                OctaTreeNode* currentNode = node;
                for (int i = 0; i < 8; ++i) {
                    Box newBox = box;
                    newBox.min = box.min + (box.max - box.min) * (Vec3f(i % 2, i / 2, i / 4)) / 2;
                    std::vector<int> v_ids_new_new;
                    for (auto id: v_ids_new) {
                        if (newBox.contains(scene.vertex_data[id-1])) {
                            v_ids_new_new.push_back(id);
                        }
                    }
                    currentNode->nextNode = build(v_ids_new_new, scene, newBox, depth + 1);
                    currentNode = currentNode->nextNode;
                }
            }



            return node;
        }
        return nullptr;
    }
};

class OctaTree {
public:
    Mesh &mesh;
    Scene &scene;
    OctaTreeNode *root;

    OctaTree(Mesh &mesh, Scene &scene) : mesh(mesh), scene(scene) {
        std::vector<int> v_ids;
        for (auto face: mesh.faces) {
            v_ids.push_back(face.v0_id);
            v_ids.push_back(face.v1_id);
            v_ids.push_back(face.v2_id);
        }
        root = OctaTreeNode::build(v_ids, scene, scene.getBoundingBox(v_ids), 0);
    }
};

class RayTracer {
public:
    RayTracer(parser::Scene &scene) {
        this->scene = scene;
    }


    void rayTrace() {
        std::vector<MySphere> meshBoundingSpheres;
        std::vector<OctaTree> meshTrees;
        for (auto &mesh: scene.meshes) {
            meshBoundingSpheres.push_back(scene.getBoundingSphere(mesh));
            meshTrees.push_back(OctaTree(mesh, scene));
        }
        int boundingMeshHit = 0;
        int boundingMeshMiss = 0;

        for (auto camera: scene.cameras) {
            auto *image = new unsigned char[camera.image_width * camera.image_height * 3];
            int imagePtr = 0;
            for (int j = 0; j < camera.image_height; j++) {
                std::cout << 100 * j / (double) camera.image_height << '%' << std::endl;
                for (int i = 0; i < camera.image_width; i++) {
                    Ray eyeRay = generateEyeRay(camera, i, j);

                    auto raytracedColor = scene.background_color;

                    if (j == 350 && i == 350) {
                        int a = 45;
                    }
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
                        if (eyeRay.intersects(scene, meshBoundingSpheres[meshIndex]).exists) {
                            boundingMeshHit++;
                            auto &mesh = scene.meshes[meshIndex];
                            for (auto &face: mesh.faces) {
                                auto intersectionPoint = eyeRay.intersects(scene, face, mesh.material_id);
                                if (intersectionPoint.exists) {

                                    raytracedColor.x = 0;
                                    raytracedColor.y = 0;
                                    raytracedColor.z = 255;


                                }
                            }
                        } else {
                            boundingMeshMiss++;
                        }

                    }

                    image[imagePtr++] = raytracedColor.x;
                    image[imagePtr++] = raytracedColor.y;
                    image[imagePtr++] = raytracedColor.z;
                }
            }
            image[3] = 100;
            image[4] = 100;
            image[5] = 100;
            write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);
            std::cout << camera.image_name << std::endl;
            printf("meshBoundingSphereHits %%%.1f \n", 100 * boundingMeshHit / (double) boundingMeshMiss);

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

    auto begin = std::chrono::high_resolution_clock::now();
    RayTracer raytracer(scene);
    raytracer.rayTrace();
    auto end = std::chrono::high_resolution_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("\nRendered in %.3f seconds.\n", elapsed.count() * 1e-9);
}
