#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <list>

namespace parser {

    typedef unsigned char Pixel[3];
    typedef Pixel *Image;


    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f {
        float x, y, z;


        Vec3f operator+(const Vec3f &v) const {
            return Vec3f{x + v.x, y + v.y, z + v.z};
        }

        Vec3f operator*(float f) const {
            return Vec3f{x * f, y * f, z * f};
        }

        float operator*(const Vec3f &v) const {
            return x * v.x + y * v.y + z * v.z;
        }

        Vec3f operator-(const Vec3f &v) const {
            return Vec3f{x - v.x, y - v.y, z - v.z};
        }

        Vec3f operator-() const {
            return Vec3f{-x, -y, -z};
        }

        Vec3f crossProduct(const Vec3f &v) const {
            return Vec3f{y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
        }

        Vec3f dotWithoutSum(const Vec3f &v) const {
            return Vec3f{x * v.x, y * v.y, z * v.z};
        }

        float distance(const Vec3f &v) const {
            return sqrt(pow(x - v.x, 2) + pow(y - v.y, 2) + pow(z - v.z, 2));
        }

        // subscript operator
        float &operator[](int i) {
            switch (i) {
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                default:
                    return x;
            }
        }

        Vec3f operator/(float f) const {
            return Vec3f{x / f, y / f, z / f};
        }

        Vec3f normalize() const {
            float vLength = length();
            return Vec3f{x / vLength, y / vLength, z / vLength};
        }

        float length() const {
            return sqrt(x * x + y * y + z * z);
        }

        void clamp(float a, float b) {
            x = std::max(a, std::min(x, b));
            y = std::max(a, std::min(y, b));
            z = std::max(a, std::min(z, b));
        }


    };

    struct Vec3i {
        int x, y, z;

        void toPixel(Pixel &pixel) {
            pixel[0] = x;
            pixel[1] = y;
            pixel[2] = z;
        }

        int &operator[](int i) {
            switch (i) {
                case 0:
                    return x;
                case 1:
                    return y;
                case 2:
                    return z;
                default:
                    return x;
            }
        }

        Vec3i operator+(const Vec3f &v) const {
            return Vec3i{(int) (x + v.x), (int) (y + v.y), (int) (z + v.z)};
        }

        Vec3i operator+(const Vec3i &v) const {
            return Vec3i{x + v.x, y + v.y, z + v.z};
        }

        Vec3i& operator+=(const Vec3i& rhs){
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        Vec3i& operator+=(const Vec3f& rhs){
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        Vec3i dotWithoutSum(const Vec3f &v) const {
            return Vec3i{(int) (x * v.x), (int) (y * v.y), (int) (z * v.z)};
        }

        void clamp(int a, int b) {
            x = std::max(a, std::min(x, b));
            y = std::max(a, std::min(y, b));
            z = std::max(a, std::min(z, b));
        }

        bool allGreaterEqualTo(int a) {
            return x >= a && y >= a && z >= a;
        }
    };

    struct Vec4f {
        float x, y, z, w;
    };

    struct Camera {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material {
        bool is_mirror;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Sphere {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct MySphere {
        Vec3f c;
        float r;

        MySphere() = default;

        MySphere(Vec3f c, float r) : c(c), r(r) {}
    };

    struct Box {
        Vec3f min, max;

        Box() = default;

        Box(Vec3f min, Vec3f max) : min(min), max(max) {}

        bool contains(Vec3f point) {
            return point.x >= min.x && point.x <= max.x &&
                   point.y >= min.y && point.y <= max.y &&
                   point.z >= min.z && point.z <= max.z;
        }
    };

    struct Mesh {
        int material_id;
        std::vector<Face> faces;
    };

    struct Triangle {
        int material_id;
        Face indices;

        Triangle() = default;

        Triangle(int material_id, Face indices) : material_id(material_id), indices(indices) {}
    };


    struct Scene {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        //Functions
        void loadFromXml(const std::string &filepath);


        Box getBoundingBox(std::list<Triangle> &faces) {
            Vec3f min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                         std::numeric_limits<float>::max()};
            Vec3f max = {-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
                         -std::numeric_limits<float>::max()};
            for (auto &face: faces) {
                for (auto &vertexId: {face.indices.v0_id, face.indices.v1_id, face.indices.v2_id}) {
                    auto &vertex = vertex_data[vertexId - 1];
                    if (vertex.x < min.x) {
                        min.x = vertex.x;
                    }
                    if (vertex.y < min.y) {
                        min.y = vertex.y;
                    }
                    if (vertex.z < min.z) {
                        min.z = vertex.z;
                    }
                    if (vertex.x > max.x) {
                        max.x = vertex.x;
                    }
                    if (vertex.y > max.y) {
                        max.y = vertex.y;
                    }
                    if (vertex.z > max.z) {
                        max.z = vertex.z;
                    }
                }
            }
            return {min, max};
        }

        void extendBoundingBox(Box &box, std::list<Sphere> &sphereList) {
            for (auto &sphere: sphereList) {
                auto &center = vertex_data[sphere.center_vertex_id - 1];
                for (int axis = 0; axis < 3; ++axis) {
                    if (center[axis] - sphere.radius < box.min[axis]) {
                        box.min[axis] = center[axis] - sphere.radius;
                    }
                    if (center[axis] + sphere.radius > box.max[axis]) {
                        box.max[axis] = center[axis] + sphere.radius;
                    }

                }

            }
        }

        Vec3f getCenter(Face &face) {
            return (vertex_data[face.v0_id - 1] + vertex_data[face.v1_id - 1] + vertex_data[face.v2_id - 1]) / 3;

        }
    };
}

#endif
