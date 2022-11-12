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
#include <mutex>
#include <condition_variable>

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

struct Intersection {
    float tSmall, tLarge;
    Vec3f normal;
    int material_id;
    bool exists;
};

struct BoxIntersection {
    float t;
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


    Intersection intersects(Sphere sphere) {
        Intersection result = {-1, -1, {0, 0, 0}, sphere.material_id, false};
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

    // source: book. also taking advantage of floating point arithmetic,
    // when a direction component is 0, result is intinity but it does not affect the result. thus we don't need branches
    // also division is slower than multiplication, so we cache one over direction
    BoxIntersection intersects(Box box) {
        auto tx1 = (box.min.x - origin.x) * oneOverDirection.x;
        auto tx2 = (box.max.x - origin.x) * oneOverDirection.x;

        auto tmin = std::min(tx1, tx2);
        auto tmax = std::max(tx1, tx2);

        auto ty1 = (box.min.y - origin.y) * oneOverDirection.y;
        auto ty2 = (box.max.y - origin.y) * oneOverDirection.y;

        tmin = std::max(tmin, std::min(ty1, ty2));
        tmax = std::min(tmax, std::max(ty1, ty2));

        auto tz1 = (box.min.z - origin.z) * oneOverDirection.z;
        auto tz2 = (box.max.z - origin.z) * oneOverDirection.z;

        tmin = std::max(tmin, std::min(tz1, tz2));
        tmax = std::min(tmax, std::max(tz1, tz2));

        if (tmax >= std::max(0.0f, tmin)) {
            return {tmin, true};
        } else {
            return {-1, false};
        }

    }


    Intersection intersects(Scene &scene, Triangle &triangle) {
        Intersection point = {-1, -1, triangle.normal, triangle.material_id, false};
        auto &a = scene.vertex_data[triangle.indices.v0_id - 1];
        auto &b = scene.vertex_data[triangle.indices.v1_id - 1];
        auto &c = scene.vertex_data[triangle.indices.v2_id - 1];

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

    Intersection getFirstIntersection(Scene &scene, BVHTree &tree) {
        Intersection firstIntersection = {-1, -1, {-1, -1}, -1, false};
        auto stack = startTraversing();
        float tMax = std::numeric_limits<float>::max();


        while (!stack.empty()) {
            auto node = stack.top();
            stack.pop();
            auto intersection = intersects(tree.nodes[node].box);

            if (intersection.exists && intersection.t <= tMax) {
                if (!tree.nodes[node].isLeaf()) {
                    if (direction[tree.nodes[node].axis] > 0) {
                        stack.push(tree.nodes[node].rightIndex);
                        stack.push(node+1);
                    } else {
                        stack.push(node+1);
                        stack.push(tree.nodes[node].rightIndex);
                    }

                } else {
                    for (auto triangle: tree.nodes[node].triangles) {
                        auto intersectionPoint = intersects(scene, triangle);
                        if (intersectionPoint.exists) {
                            if (intersectionPoint.tSmall < firstIntersection.tSmall || firstIntersection.tSmall == -1) {
                                firstIntersection = intersectionPoint;
                                tMax = firstIntersection.tSmall;
                            }
                        }
                    }
                    for (auto &sphere: tree.nodes[node].spheres) {
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

    Intersection getAnyIntersectionUntilT(Scene &scene, BVHTree &tree, float t) {
        Intersection firstIntersection = {-1, -1, {0}, -1, false};
        auto stack = startTraversing();

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

    std::stack<int> startTraversing() {
        std::stack<int> stack;
        if (tree.root) {
            stack.push(0);
        }
        return stack;
    }


    BVHNode *traverse(std::stack<int> &stack) {
        auto node = stack.top();
        stack.pop();
        if (intersects(tree.nodes[node].box).exists) {
            if (!tree.nodes[node].isLeaf()) {
                if (direction[tree.nodes[node].axis] > 0) {
                    stack.push(tree.nodes[node].rightIndex);
                    stack.push(node+1);
                } else {
                    stack.push(node+1);
                    stack.push(tree.nodes[node].rightIndex);
                }
            }
            return &(tree.nodes[node]);
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
    Camera *currentCamera;
    Image currentImage;
    EyeRayGenerator eyeRayGenerator;

    explicit RayTracer(parser::Scene &scene) : scene(scene), tree(scene), eyeRayGenerator(scene, tree) {
        std::vector<Triangle> triangles(scene.triangles.begin(), scene.triangles.end());
        for (auto &mesh: scene.meshes) {
            for (auto &face: mesh.faces) {
                triangles.emplace_back(mesh.material_id, face);
            }
        }
        for (auto &triangle: triangles) {
            auto &a = scene.vertex_data[triangle.indices.v0_id - 1];
            auto &b = scene.vertex_data[triangle.indices.v1_id - 1];
            auto &c = scene.vertex_data[triangle.indices.v2_id - 1];
            triangle.normal = ((b - a).crossProduct(c - a)).normalize();
            triangle.center = (a + b + c) / 3;
        }
        tree.build(triangles);
    }

    void renderWithMultipleThreads(int threadNumber, int totalThreads, std::mutex *m, std::condition_variable *cv, bool *ready, int *lastRenderedRow, std::vector<bool> *renderedRows) {
        for (int rowNum = threadNumber; rowNum < currentCamera->image_height; rowNum+=totalThreads){
            for (int colNum = 0; colNum < currentCamera->image_width; colNum++) {
                Ray eyeRay = eyeRayGenerator.generate(rowNum, colNum);
                auto raytracedColor = rayTrace(eyeRay);
                raytracedColor.toPixel(currentImage[currentCamera->image_width * rowNum + colNum]);
            }
            (*renderedRows)[rowNum] = true;
            if (rowNum % 64 == 63) {
                *lastRenderedRow = rowNum;
                std::lock_guard<std::mutex> lk(*m);
/*
                std::cout << "Thread " << threadNumber << " finished row " << rowNum << std::endl;
*/
                *ready = true;
                cv->notify_all();
            }

        }
    }

    void render(Camera &camera) {
        auto image = new Pixel[camera.image_width * camera.image_height];
        currentCamera = &camera;
        currentImage = image;
        eyeRayGenerator.init(currentCamera);
        std::mutex m;
        std::condition_variable cv;
        int lastRenderedRow = 0;
        bool ready = false;
        std::vector<bool> renderedRows(camera.image_height, false);


        auto processor_count = std::thread::hardware_concurrency();
        if (processor_count == 0) {
            processor_count = 8; // todo: try out different thread numbers for inek. maybe inek does not have 4 cores
        }
        processor_count--;

        std::vector<std::thread> threads;
        threads.reserve(processor_count);
        std::cout << "Rendering " << camera.image_name << " with " << processor_count << " threads (cores) and writing to disk with 1 thread" << std::endl;

        std::thread ppmWriterThread([&](){
            PPMWriter writer(camera.image_name.c_str(), (unsigned char *) image, camera.image_width, camera.image_height);
            while (writer.lastWrittenRow < writer.height-1) {
                std::unique_lock<std::mutex> lk(m);
                cv.wait(lk, [&](){return ready;});
/*
                std::cout << "Writing row from " << writer.lastWrittenRow << "to " << lastRenderedRow << std::endl;
*/
                int lastRenderedRowCopy = lastRenderedRow;
                ready = false;
                lk.unlock();
                for (; writer.lastWrittenRow < lastRenderedRowCopy;) {
                    if (renderedRows[writer.lastWrittenRow+1]) {
                        writer.writeNextRow();
                    } else {
                        break;
                    }
                }
                fflush(writer.outfile);

            }

        });



        auto begin2 = std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < processor_count; i++) {
            threads.emplace_back(&RayTracer::renderWithMultipleThreads, this, i, processor_count, &m, &cv, &ready, &lastRenderedRow, &renderedRows);
        }
        for (auto &thread: threads) {
            thread.join();
        }
        {
            std::lock_guard<std::mutex> lk(m);
            ready = true;
            lastRenderedRow = camera.image_height-1;
            cv.notify_all();
        }

        ppmWriterThread.join();

    }

    Vec3f rayTrace(Ray &ray, int depth = 0) {
        Vec3f rayTracedColor = {0, 0, 0};
        if (depth > ray.scene.max_recursion_depth) {
            return rayTracedColor;
        }
        Intersection intersection = ray.getFirstIntersection(scene, tree);
        if (intersection.exists) {
            auto &material = scene.materials[intersection.material_id - 1];

            auto ambient = material.ambient.dotWithoutSum(scene.ambient_light);
            rayTracedColor += ambient;

            Vec3f intersectionPnt = ray.getPoint(intersection.tSmall) + intersection.normal * scene.shadow_ray_epsilon;

            for (auto &light: scene.point_lights) {
                auto lightDistance = (light.position - intersectionPnt).length();
                auto lightRayDirection = (light.position - intersectionPnt).normalize();
                auto lightRayDirectionReal = (light.position - ray.getPoint(intersection.tSmall)).normalize();
                auto lightRay = Ray(intersectionPnt, lightRayDirection, scene, tree);
                auto lightRayIntersection = lightRay.getAnyIntersectionUntilT(scene, tree, lightDistance);

                if (!lightRayIntersection.exists) {
                    float cosTheta = lightRayDirectionReal * intersection.normal;
                    auto receivedIrradiance = light.intensity / (lightDistance * lightDistance);


                    float theta = acos(cosTheta) * 180 / 3.1415;
                    if (theta <= 90.01) { // without lightRayDirectionReal and 0.01, this is false instead of true
                        auto h = (lightRay.direction + -ray.direction.normalize()).normalize();
                        float cosaToTheP = pow(std::max(0.0f, intersection.normal.normalize() * h),
                                               material.phong_exponent);
                        auto specular = (material.specular * cosaToTheP).dotWithoutSum(receivedIrradiance);
                        rayTracedColor += specular;
                    }

                    auto diffuse = (material.diffuse * clampFloat(cosTheta, 0, 1)).dotWithoutSum(
                            receivedIrradiance);
                    rayTracedColor += diffuse;


                }

            }


            if (material.is_mirror) {
                ray.direction = ray.direction.normalize();
                intersection.normal = intersection.normal.normalize();
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
                return {(float) scene.background_color.x, (float) scene.background_color.y,
                        (float) scene.background_color.z};
            }

        }
        return rayTracedColor.clamp(0, std::numeric_limits<float>::max());
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
            rayTracer.render(camera);
        }
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2);
    printf("Rendered in %.3f seconds.\n", elapsed2.count() * 1e-9 / renderCount);
    printf("Total: %.3f seconds.\n", (elapsed2.count() * 1e-9 / renderCount) + (elapsed1.count() * 1e-9));
}
