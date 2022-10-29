#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <chrono>


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

class RayTracer {
public:
    RayTracer(parser::Scene &scene) {
        this->scene = scene;
    }


    void rayTrace() {
        std::vector<MySphere> meshBoundingSpheres;
        for (auto &mesh : scene.meshes) {
            meshBoundingSpheres.push_back(scene.getBoundingSphere(mesh));
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

                    for (auto &sphere: scene.spheres) {
                        auto intersectionPoint = eyeRay.intersects(scene, { scene.vertex_data[sphere.center_vertex_id - 1], sphere.radius});
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
                            auto & mesh = scene.meshes[meshIndex];
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
            write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);
            std::cout << camera.image_name << std::endl;
            printf("meshBoundingSphereHits %%%.1f \n", 100 * boundingMeshHit/(double)boundingMeshMiss);

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
