#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>

typedef unsigned char RGB[3];
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


    IntersectionPoint intersects(Scene &scene, Sphere &sphere) {
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
            return {t1, t2, true};
        }
        return {0, 0, false};

    }

    bool intersects(Scene &scene, Face &face, int material_id) {
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

        return alpha >= 0 &&
               beta >= 0 &&
               gamma >= 0 &&
               t >= tMin;


    }
};

class RayTracer {
public:
    RayTracer(parser::Scene &scene) {
        this->scene = scene;
    }


    void rayTrace() {
        for (auto camera: scene.cameras) {
            unsigned char *image = new unsigned char[camera.image_width * camera.image_height * 3];
            int imagePtr = 0;
            for (int j = 0; j < camera.image_height; j++) {
                std::cout << 100 * j / (double) camera.image_height << '%' << std::endl;
                for (int i = 0; i < camera.image_width; i++) {


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

                    Ray eyeRay(e, s - e);

                    RGB raytracedColor = {0, 0, 0};

                    for (auto sphere: scene.spheres) {
                        auto intersectionPoint = eyeRay.intersects(scene, sphere);
                        if (intersectionPoint.exists) {
                            raytracedColor[0] = 255;
                            raytracedColor[1] = 0;
                            raytracedColor[2] = 0;
                        }
                    }
                    for (auto triangle: scene.triangles) {
                        if (eyeRay.intersects(scene, triangle.indices, triangle.material_id)) {
                            raytracedColor[0] = 255;
                            raytracedColor[1] = 255;
                            raytracedColor[2] = 0;
                        }
                    }
                    for (auto mesh: scene.meshes) {
                        for (auto face: mesh.faces) {
                            if (eyeRay.intersects(scene, face, mesh.material_id)) {
                                raytracedColor[0] = 0;
                                raytracedColor[1] = 0;
                                raytracedColor[2] = 255;
                            }
                        }
                    }
                    // todo: raytracing

                    image[imagePtr++] = raytracedColor[0];
                    image[imagePtr++] = raytracedColor[1];
                    image[imagePtr++] = raytracedColor[2];
                }
            }
            write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);
            std::cout << camera.image_name;
        }

    }


private:
    parser::Scene scene;

};

int main(int argc, char *argv[]) {
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

    RayTracer raytracer(scene);
    raytracer.rayTrace();
}
