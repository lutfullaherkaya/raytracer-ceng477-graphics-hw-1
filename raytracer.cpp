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
    float tSmall, tLarge;
    Vec3f normal;
    int material_id;
    bool exists;
};

class Ray {
public:
    Vec3f getPoint(float t) { // t >= 0
        return origin + direction * t;
    }

    Vec3f origin;
    Vec3f direction;


    Ray(Vec3f origin, Vec3f direction) : origin(origin), direction(direction) {
        direction = direction.normalize();
    }


    IntersectionPoint intersects(Scene &scene, Sphere sphere) {
        IntersectionPoint result = {-1, -1, {0,0,0}, sphere.material_id, false};
        auto c = scene.vertex_data[sphere.center_vertex_id - 1];
        auto r = sphere.radius;
        auto d = direction;
        auto o = origin;

        auto B = 2 * (d * (o - c));
        auto A = d * d;
        auto C = (o - c) * (o - c) - r * r;
        auto discriminant = B * B - 4 * A * C;

        if (discriminant >= 0) {
            /*float t1 = (-d * (o - c) - sqrt(discriminant)) / (2 * (d * d));*/
            float t1 = (-B - sqrt(discriminant)) / (2 * A);
            float t2 = (-B + sqrt(discriminant)) / (2 * A);
            if (t1 < 0 && t2 < 0) {
                return result;
            }
            result.exists = true;
            result.tSmall = t1;
            result.tLarge = t2;
            result.normal = ((getPoint(t1) - c)/r).normalize();
            return result;
        }
        return result;

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
            return {-1, -1, false};

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
        if (tmin < 0 && tmax < 0) {
            return {-1, -1, {0,0,0}, -1, false};
        }
        return {tmin, tmax, {0,0,0}, -1, true};

    }

    IntersectionPoint intersects(Scene &scene, Triangle &triangle) {
        IntersectionPoint point = {-1, -1, {-1,-1,-1}, triangle.material_id, false};
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
            point.tSmall = t;
        }
        return point;

    }

    IntersectionPoint getFirstIntersection(Scene &scene, BVHTree &tree) {
        IntersectionPoint firstIntersection = {-1, -1, -1, false};
        std::stack<BVHNode *> stack;
        stack.push(tree.root);
        while (!stack.empty()) {
            auto node = stack.top();
            stack.pop();
            if (intersects(scene, node->box).exists) {
                if (node->left) {
                    stack.push(node->left);
                }
                if (node->right) {
                    stack.push(node->right);
                }
                if (!node->left && !node->right) {
                    for (auto triangle: node->faces) {
                        auto intersectionPoint = intersects(scene, triangle);
                        if (intersectionPoint.exists) {
                            if (intersectionPoint.tSmall < firstIntersection.tSmall || firstIntersection.tSmall == -1) {
                                firstIntersection = intersectionPoint;
                            }
                        }
                    }
                    for (auto &sphere: node->spheres) {
                        auto intersectionPoint = intersects(scene, sphere);
                        if (intersectionPoint.exists) {
                            if (intersectionPoint.tSmall < firstIntersection.tSmall || firstIntersection.tSmall == -1) {
                                firstIntersection = intersectionPoint;
                            }
                        }
                    }
                }

            }
        }
        return firstIntersection;
    }
};


class RayTracer {
public:
    parser::Scene scene;
    BVHTree tree;

    explicit RayTracer(parser::Scene &scene) : scene(scene), tree(scene) {
        std::list<Triangle> triangles(scene.triangles.begin(), scene.triangles.end());
        for (auto &mesh: scene.meshes) {
            for (auto &face: mesh.faces) {
                triangles.emplace_back(mesh.material_id, face);
            }
        }
        tree.build(triangles);
    }


    Image render(Camera &camera) {
        auto image = new Pixel[camera.image_width * camera.image_height];
        for (int rowNum = 0, imagePtr = 0; rowNum < camera.image_height; rowNum++) {
            for (int colNum = 0; colNum < camera.image_width; colNum++) {
                auto raytracedColor = rayTrace(camera, rowNum, colNum);
                raytracedColor.toPixel(image[imagePtr++]);
            }
        }
        return image;

    }

    Vec3i rayTrace(Camera &camera, int rowNum, int colNum) {
        Ray eyeRay = generateEyeRay(camera, rowNum, colNum);
        auto raytracedColor = scene.background_color;

        IntersectionPoint firstIntersection = eyeRay.getFirstIntersection(scene, tree);
        if (firstIntersection.exists) {
            auto &material = scene.materials[firstIntersection.material_id-1];
            Vec3f diffuse = {0, 0, 0};
            for (auto &light: scene.point_lights) {
                float cosTheta = 0;
                auto lightRayDirection = (light.position - eyeRay.getPoint(firstIntersection.tSmall)).normalize();
                auto lightDistance = (light.position - eyeRay.getPoint(firstIntersection.tSmall)).length();
                auto lightRay = Ray(eyeRay.getPoint(firstIntersection.tSmall), lightRayDirection);
                auto lightIntersection = lightRay.getFirstIntersection(scene, tree);

                bool lightIsObstructed = lightIntersection.exists && (lightRay.direction * lightIntersection.tSmall).length() < lightDistance;
                if (!lightIsObstructed) {
                    // todo: hatayi buldum. kurenin ustunden isin yollayinca kure kendi kendine engel oluyor floating point seyinden oturu.
                    cosTheta = lightRayDirection * firstIntersection.normal.normalize();
                    if (cosTheta < 0) {
                        cosTheta = 0;
                    }
                    if (cosTheta > 1) {
                        cosTheta = 1;
                    }

                    std::cout << cosTheta;
                    for (int axis = 0; axis < 3; ++axis) {
                        diffuse[axis] += material.diffuse[axis] * cosTheta * light.intensity[axis] / (lightDistance * lightDistance); 
                    }
                    float ambient = material.ambient[0] * scene.ambient_light[0];

                }
            }
            for (int axis = 0; axis < 3; ++axis) {
                float ambient = material.ambient[axis] * scene.ambient_light[axis];
                raytracedColor[axis] = diffuse[axis] + ambient;
                if (raytracedColor[axis] > 255) {
                    raytracedColor[axis] = 255;
                }


            }
            
        }
        return raytracedColor;
    }

    static Ray generateEyeRay(Camera &camera, int rowNum, int colNum) {
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

        float su = (colNum + 0.5) * (r - l) / nx;
        float sv = (rowNum + 0.5) * (t - b) / ny;

        auto s = q + u * su - v * sv;

        return {e, s - e};
    }


};

int main(int argc, char *argv[]) {
    parser::Scene scene;
    scene.loadFromXml(argv[1]);


    auto begin1 = std::chrono::high_resolution_clock::now();
    RayTracer rayTracer(scene);
    auto end1 = std::chrono::high_resolution_clock::now();
    auto elapsed1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1);
    printf("Planted trees in %.3f seconds.\n", elapsed1.count() * 1e-9);


    auto begin2 = std::chrono::high_resolution_clock::now();
    int renderCount = 1; // todo: make it 1. 10 is for performance measurement
    for (int i = 0; i < renderCount; ++i) {
        for (auto camera: scene.cameras) {
            auto image = rayTracer.render(camera);
            write_ppm(camera.image_name.c_str(), (unsigned char *)image, camera.image_width, camera.image_height);
        }
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2);
    printf("Rendered in %.3f seconds.\n", elapsed2.count() * 1e-9 / renderCount);
}
