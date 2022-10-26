#include <iostream>
#include "parser.h"
#include "ppm.h"

typedef unsigned char RGB[3];

class RayTracer
{
public:
    RayTracer(parser::Scene &scene) {
        this->scene = scene;
    }
    void rayTrace()
    {
        for (auto camera : scene.cameras)
        {
            unsigned char *image = new unsigned char[camera.image_width * camera.image_height * 3];
            for (int i = 0; i < camera.image_height; i++)
            {
                for (int j = 0; j < camera.image_width; j++)
                {
                    RGB raytracedColor = {0, 0, 0};

                    // do raytrace

                    image[i++] = raytracedColor[0];
                    image[i++] = raytracedColor[1];
                    image[i++] = raytracedColor[2];
                }
            }
            write_ppm(camera.image_name.c_str(), image, camera.image_width , camera.image_height);
            std::cout << camera.image_name;
        }
        
    }

private:
    parser::Scene scene;

};

int main(int argc, char *argv[])
{
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
