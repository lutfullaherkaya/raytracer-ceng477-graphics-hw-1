#include "parser.h"
#include "ppm.h"
#include "bvh.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <list>
#include <stack>
#include <utility>

using namespace parser;

float det(float m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}


struct IntersectionPoint {
    float t1, t2;
    //Vec3f normal;
    int material_id;
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


    IntersectionPoint intersects(Scene &scene, Sphere sphere) {
        auto c = scene.vertex_data[sphere.center_vertex_id - 1];
        auto r = sphere.radius;
        auto d = direction;
        auto o = origin;

        auto B = 2 * (d * (o - c));
        auto A = d * d;
        auto C = (o - c) * (o - c) - r * r;
        auto discriminant = B * B - 4 * A * C;

        if (discriminant >= 0) {
            float t1 = (-d * (o - c) - sqrt(discriminant)) / (2 * (d * d));
            float t2 = (-d * (o - c) + sqrt(discriminant)) / (2 * (d * d));
            return {t1, t2, sphere.material_id, true};
        }
        return {-1, -1, sphere.material_id, false};

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

        return {tmin, tmax, -1, true};

    }

    IntersectionPoint intersects(Scene &scene, Triangle &triangle) {
        IntersectionPoint point = {-1, -1, triangle.material_id, false};
        auto a = scene.vertex_data[triangle.indices.v0_id - 1];
        auto b = scene.vertex_data[triangle.indices.v1_id - 1];
        auto c = scene.vertex_data[triangle.indices.v2_id - 1];

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


class RayTracer {
public:
    RayTracer(parser::Scene &scene, std::list<Triangle> triangles) : triangleTree(std::move(triangles), scene),
                                                                     scene(scene) {}


    BVHTree triangleTree;

    Image rayTrace(Camera &camera) {
        auto image = new unsigned char[camera.image_width * camera.image_height * 3];
        int imagePtr = 0;
        for (int j = 0; j < camera.image_height; j++) {
            std::cout << 100 * j / (double) camera.image_height << '%' << std::endl;
            for (int i = 0; i < camera.image_width; i++) {
                Ray eyeRay = generateEyeRay(camera, i, j);

                auto raytracedColor = scene.background_color;
                // checks for both triangles and meshes
                std::stack<BVHNode *> stack;
                stack.push(triangleTree.root);
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
                            for (auto triangle: node->faces) {
                                auto intersectP = eyeRay.intersects(scene, triangle);
                                if (intersectP.exists) {
                                    raytracedColor.x = 0;
                                    raytracedColor.y = 0;
                                    raytracedColor.z = 255;
                                }
                            }
                            for (auto &sphere: node->spheres) {
                                auto intersectionPoint = eyeRay.intersects(scene, sphere);
                                if (intersectionPoint.exists) {
                                    raytracedColor.x = 255;
                                    raytracedColor.y = 0;
                                    raytracedColor.z = 0;
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
        return image;

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


    auto begin1 = std::chrono::high_resolution_clock::now();
    std::list<Triangle> triangles(scene.triangles.begin(), scene.triangles.end());
    for (auto &mesh: scene.meshes) {
        for (auto &face: mesh.faces) {
            triangles.emplace_back(mesh.material_id, face);
        }
    }
    RayTracer rayTracer(scene, triangles);
    auto end1 = std::chrono::high_resolution_clock::now();


    auto begin2 = std::chrono::high_resolution_clock::now();
    for (auto camera: scene.cameras) {
        auto image = rayTracer.rayTrace(camera);
        write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);
        std::cout << camera.image_name << std::endl;
    }
    auto end2 = std::chrono::high_resolution_clock::now();

    auto elapsed1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1);
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2);
    printf("\nPlanted trees in %.3f seconds.\n", elapsed1.count() * 1e-9);
    printf("Rendered in %.3f seconds.\n", elapsed2.count() * 1e-9);
}
