#include "parser.h"
#include "ppm.h"
#include "bvh.h"

#include <iostream>
#include <cmath>
#include <chrono>
#include <list>
#include <stack>
#include <utility>
#include <thread>

using namespace parser;

float det(float m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

float clampFloat(float x, float min, float max) {
    return std::max(min, std::min(max, x));
}


#define T_MIN_EPSILON 0.001f

struct IntersectionPoint {
    float tSmall, tLarge;
    Vec3f normal;
    int material_id;
    bool exists;
};

bool floatEquals(float a, float b) {
    return fabs(a - b) < 0.001f;
}

class Ray {
public:
    Vec3f getPoint(float t) { // t >= 0
        return origin + direction * t;
    }

    const Vec3f origin;
    Vec3f direction;
    const Vec3f oneOverDirection; // multiplication is faster than division so we cache this.
    Scene &scene;
    BVHTree &tree;


    Ray(Vec3f origin, Vec3f direction, Scene &scene, BVHTree &tree) : origin(origin), direction(direction),
                                                                      oneOverDirection{1 / direction.x, 1 / direction.y,
                                                                                       1 / direction.z},
                                                                      scene(scene), tree(tree) {
        direction = direction.normalize();

    }


    IntersectionPoint intersects(Sphere sphere) {
        IntersectionPoint result = {-1, -1, {0, 0, 0}, sphere.material_id, false};
        auto c = scene.vertex_data[sphere.center_vertex_id - 1];
        auto r = sphere.radius;
        auto d = direction;
        auto o = origin;

        auto B = 2 * (d * (o - c));
        auto A = d * d;
        auto C = (o - c) * (o - c) - r * r;
        auto discriminant = B * B - 4 * A * C;

        if (discriminant >= 0) {
            float t1 = (-B - sqrt(discriminant)) / (2 * A);
            float t2 = (-B + sqrt(discriminant)) / (2 * A);
            if (t1 < 0 && t2 < 0) {
                return result;
            }
            result.exists = true;
            result.tSmall = t1;
            result.tLarge = t2;
            result.normal = ((getPoint(t1) - c) / r).normalize();
            return result;
        }
        return result;

    }



    // source: book
    IntersectionPoint intersects(Box box) {
        float tmin = (box.min.x - origin.x) / direction.x;
        float tmax = (box.max.x - origin.x) / direction.x;

        if (tmin > tmax) std::swap(tmin, tmax);

        float tymin = (box.min.y - origin.y) / direction.y;
        float tymax = (box.max.y - origin.y) / direction.y;

        if (tymin > tymax) std::swap(tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax))
            return {-1, -1, {0}, -1, false};

        if (tymin > tmin)
            tmin = tymin;

        if (tymax < tmax)
            tmax = tymax;

        float tzmin = (box.min.z - origin.z) / direction.z;
        float tzmax = (box.max.z - origin.z) / direction.z;

        if (tzmin > tzmax) std::swap(tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax))
            return {0, 0, {0}, -1, false};

        if (tzmin > tmin)
            tmin = tzmin;

        if (tzmax < tmax)
            tmax = tzmax;
        if (tmin < 0 && tmax < 0) {
            return {-1, -1, {0, 0, 0}, -1, false};
        }
        return {tmin, tmax, {0, 0, 0}, -1, true};

    }


    IntersectionPoint intersects(Scene &scene, Triangle &triangle) {
        IntersectionPoint point = {-1, -1, {-1, -1, -1}, triangle.material_id, false};
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
            point.normal = ((b - a).crossProduct(c - a)).normalize();
        }
        return point;

    }

    IntersectionPoint getFirstIntersection(Scene &scene, BVHTree &tree) {
        IntersectionPoint firstIntersection = {-1, -1, {-1, -1}, -1, false};
        std::stack<BVHNode *> stack = startTraversing();
        float tMax = std::numeric_limits<float>::max();


        while (!stack.empty()) {
            auto node = stack.top();
            stack.pop();
            auto intersection = intersects(node->box);

            if (intersection.exists && intersection.tSmall <= tMax) {
                if (!node->isLeaf()) {
                    if (direction[node->axis] > 0) {
                        stack.push(node->right);
                        stack.push(node->left);
                    } else {
                        stack.push(node->left);
                        stack.push(node->right);
                    }

                } else {
                    for (auto triangle: node->triangles) {
                        auto intersectionPoint = intersects(scene, triangle);
                        if (intersectionPoint.exists) {
                            if (intersectionPoint.tSmall < firstIntersection.tSmall || firstIntersection.tSmall == -1) {
                                firstIntersection = intersectionPoint;
                                tMax = firstIntersection.tSmall;
                            }
                        }
                    }
                    for (auto &sphere: node->spheres) {
                        auto intersectionPoint = intersects(sphere);
                        if (intersectionPoint.exists) {
                            if (intersectionPoint.tSmall < firstIntersection.tSmall || firstIntersection.tSmall == -1) {
                                firstIntersection = intersectionPoint;
                                tMax = firstIntersection.tSmall;
                            }
                        }
                    }

                }
            }


        }

        return firstIntersection;
    }

    IntersectionPoint getAnyIntersectionUntilT(Scene &scene, BVHTree &tree, float t) {
        IntersectionPoint firstIntersection = {-1, -1, {0}, -1, false};
        std::stack<BVHNode *> stack = startTraversing();

        while (!stack.empty()) {
            auto node = traverse(stack);
            if (node && node->isLeaf()) {
                for (auto triangle: node->triangles) {
                    auto intersectionPoint = intersects(scene, triangle);
                    if (intersectionPoint.exists) {
                        if (intersectionPoint.tSmall < t) { // todo: maybe float error for t
                            return intersectionPoint;
                        }
                    }
                }
                for (auto &sphere: node->spheres) {
                    auto intersectionPoint = intersects(sphere);
                    if (intersectionPoint.exists) {
                        if (intersectionPoint.tSmall < t) { // todo: maybe float error for t
                            return intersectionPoint;
                        }
                    }
                }
            }
        }
        return firstIntersection;
    }

    std::stack<BVHNode *> startTraversing() {
        std::stack<BVHNode *> stack;
        stack.push(tree.root);
        return stack;
    }

    BVHNode *traverseForHitNode(std::stack<BVHNode *> &stack) {
        auto node = stack.top();
        stack.pop();
        auto intersection = intersects(node->box);
        if (intersection.exists) {
            if (!node->isLeaf()) {
                if (direction[node->axis] > 0) {
                    stack.push(node->right);
                    stack.push(node->left);
                } else {
                    stack.push(node->left);
                    stack.push(node->right);
                }

            }
            return node;
        }
        return nullptr;
    }

    BVHNode *traverse(std::stack<BVHNode *> &stack) {
        auto node = stack.top();
        stack.pop();
        if (intersects(node->box).exists) {
            if (!node->isLeaf()) {
                if (direction[node->axis] > 0) {
                    stack.push(node->right);
                    stack.push(node->left);
                } else {
                    stack.push(node->left);
                    stack.push(node->right);
                }
            }
            return node;
        }
        return nullptr;
    }

};

struct EyeRayGenerator {
    Scene &scene;
    BVHTree &tree;
    Vec3f q, u, v, e;
    float suMultiplier, svMultiplier;

    EyeRayGenerator(Scene &scene, BVHTree &tree) : scene(scene), tree(tree) {}

    void init(Camera *camera) {
        e = camera->position;
        auto w = -camera->gaze;
        auto distance = camera->near_distance;

        auto l = camera->near_plane.x;
        auto r = camera->near_plane.y;
        auto b = camera->near_plane.z;
        auto t = camera->near_plane.w;

        v = camera->up;
        u = v.crossProduct(w);
        auto nx = camera->image_width;
        auto ny = camera->image_height;

        auto m = e + -w * distance;
        q = m + u * l + v * t;

        suMultiplier = (r - l) / (float) nx;
        svMultiplier = (t - b) / (float) ny;


    }

    // copy assignment


    Ray generate(int rowNum, int colNum) {
        float su = (colNum + 0.5) * suMultiplier;
        float sv = (rowNum + 0.5) * svMultiplier;
        auto s = q + u * su - v * sv;
        return {e, s - e, scene, tree};
    }
};

class RayTracer {
public:
    parser::Scene scene;
    BVHTree tree;
    BVHTree backfaceCulledTree;
    Camera *currentCamera;
    Image currentImage;
    EyeRayGenerator eyeRayGenerator;

    explicit RayTracer(parser::Scene &scene) : scene(scene), tree(scene), backfaceCulledTree(scene), eyeRayGenerator(scene, tree) {
        std::list<Triangle> triangles(scene.triangles.begin(), scene.triangles.end());
        for (auto &mesh: scene.meshes) {
            for (auto &face: mesh.faces) {
                triangles.emplace_back(mesh.material_id, face);
            }
        }
        tree.build(triangles);



    }

    void buildTheBackfaceCulledTree() {
        std::list<Triangle> triangles(scene.triangles.begin(), scene.triangles.end());
        for (auto &mesh: scene.meshes) {
            for (auto &face: mesh.faces) {
                triangles.emplace_back(mesh.material_id, face);
            }
        }


        std::list<Triangle> trianglesBackfaceCulled;
        // filter triangles by backface culling
        for (auto& triangle: triangles) {
            if (triangleFacesCamera(triangle)) {
                trianglesBackfaceCulled.push_back(triangle);
            }
        }
        backfaceCulledTree.build(trianglesBackfaceCulled);
    }

    bool triangleFacesCamera(Triangle &triangle) {
        auto a = scene.vertex_data[triangle.indices.v0_id - 1];
        auto b = scene.vertex_data[triangle.indices.v1_id - 1];
        auto c = scene.vertex_data[triangle.indices.v2_id - 1];
        auto normal = ((b - a).crossProduct(c - a)).normalize();

        auto direction = a - currentCamera->position;
        return normal * direction < 0;
    }

    void renderRowsWithModulo(int threadNumber, int totalThreads) {
        for (int rowNum = threadNumber; rowNum < currentCamera->image_height; rowNum += totalThreads) {
            for (int colNum = 0; colNum < currentCamera->image_width; colNum++) {
                Ray eyeRay = eyeRayGenerator.generate(rowNum, colNum);
                auto raytracedColor = rayTrace(eyeRay);
                raytracedColor.toPixel(currentImage[currentCamera->image_width * rowNum + colNum]);
            }
        }
    }

    Image render(Camera &camera) {
        auto image = new Pixel[camera.image_width * camera.image_height];
        currentCamera = &camera;
        buildTheBackfaceCulledTree();

        currentImage = image;
        eyeRayGenerator.init(currentCamera);
        auto processor_count = std::thread::hardware_concurrency();
        if (processor_count == 0) {
            processor_count = 8;
        }
        std::vector<std::thread> threads;
        threads.reserve(processor_count);
        std::cout << "Rendering with " << processor_count << " threads (cores)..." << std::endl;
        for (unsigned int i = 0; i < processor_count; i++) {
            threads.emplace_back(&RayTracer::renderRowsWithModulo, this, i, processor_count);
        }
        for (auto &thread: threads) {
            thread.join();
        }


        return image;

    }

    Vec3i rayTrace(Ray &ray, int depth = 0) {
        Vec3i rayTracedColor = {0, 0, 0};
        if (depth > ray.scene.max_recursion_depth) {
            return rayTracedColor;
        }
        IntersectionPoint intersection;
        if (depth == 0) {
            intersection = ray.getFirstIntersection(scene, backfaceCulledTree);
        } else {
            intersection = ray.getFirstIntersection(scene, tree);
        }

        if (intersection.exists) {
            auto &material = scene.materials[intersection.material_id - 1];

            auto ambient = material.ambient.dotWithoutSum(scene.ambient_light);
            rayTracedColor += ambient;

            Vec3f intersectionPnt = ray.getPoint(intersection.tSmall) + intersection.normal * scene.shadow_ray_epsilon;

            for (auto &light: scene.point_lights) {
                if (rayTracedColor.allGreaterEqualTo(255)) {
                    break;
                }
                auto lightDistance = (light.position - intersectionPnt).length();
                auto lightRayDirection = (light.position - intersectionPnt).normalize();
                auto lightRay = Ray(intersectionPnt, lightRayDirection, scene, tree);
                auto lightRayIntersection = lightRay.getAnyIntersectionUntilT(scene, tree, lightDistance);

                if (!lightRayIntersection.exists && !rayTracedColor.allGreaterEqualTo(255)) {
                    float cosTheta = lightRayDirection * intersection.normal;
                    auto receivedIrradiance = light.intensity / (lightDistance * lightDistance);


                    float theta = acos(cosTheta) * 180 / 3.14159265358;
                    if (0 < theta && theta < 90) {
                        auto h = (lightRay.direction + -ray.direction).normalize();
                        float cosaToTheP = pow(std::max(0.0f, intersection.normal * h), material.phong_exponent);
                        auto specular = (material.specular * cosaToTheP).dotWithoutSum(receivedIrradiance);
                        rayTracedColor += specular;
                    }

                    if (!rayTracedColor.allGreaterEqualTo(255)) {
                        auto diffuse = (material.diffuse * clampFloat(cosTheta, 0, 1)).dotWithoutSum(
                                receivedIrradiance);
                        rayTracedColor += diffuse;
                    }

                }

            }


            if (material.is_mirror && !rayTracedColor.allGreaterEqualTo(255)) {
                auto reflectionCosTheta = -ray.direction * intersection.normal;
                Ray reflectionRay(intersectionPnt,
                                  ray.direction + intersection.normal * 2 * reflectionCosTheta, scene, tree);
                auto reflectedColor = rayTrace(reflectionRay, depth + 1);

                rayTracedColor = rayTracedColor + reflectedColor.dotWithoutSum(material.mirror);

            }


        } else {
            if (depth > 0) {
                return {0, 0, 0};
            } else {
                return scene.background_color;
            }

        }
        rayTracedColor.clamp(0, 255);
        return rayTracedColor;
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


    for (auto camera: scene.cameras) { // todo: sil. bunu cache misslerin performans olcumunu etkilememesi icin yapiyorum.
        auto image = rayTracer.render(camera);
        write_ppm(camera.image_name.c_str(), (unsigned char *) image, camera.image_width, camera.image_height);
    }
    auto begin2 = std::chrono::high_resolution_clock::now();
    int renderCount = 10; // todo: make it 1. 10 is for performance measurement
    for (int i = 0; i < renderCount; ++i) {
        for (auto camera: scene.cameras) {
            auto image = rayTracer.render(camera);
            write_ppm(camera.image_name.c_str(), (unsigned char *) image, camera.image_width, camera.image_height);
        }
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2);
    printf("Rendered in %.3f seconds.\n", elapsed2.count() * 1e-9 / renderCount);
}
