#include <iostream>
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];
using namespace parser;

class Ray {
    Vec3f origin;
    Vec3f direction;

    Vec3f getPoint(float t) {
        return origin + direction * t;
    }

public:
    Ray(Vec3f origin, Vec3f direction) {
        this->origin = origin;
        this->direction = direction;
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
            for (int i = 0; i < camera.image_height; i++) {
                for (int j = 0; j < camera.image_width; j++) {
                    RGB raytracedColor = {1, 2, 3};

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

                    auto m = e +  -w * distance;
                    auto q = m + u*l + v*t;

                    float su = (i + 0.5) * (r-l) / nx;
                    float sv = (j + 0.5) * (t-b) / ny;

                    auto s = q + u*su - v*sv;

                    Ray eyeRay(e, s-e);

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
