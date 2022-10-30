#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <list>

namespace parser {

    bool floatEquals(float a, float b) {
        return fabs(a - b) < std::numeric_limits<float>::epsilon();
        // todo: not: bu epsilon en kucuk epsilon ama biraz cok kucuk gelebilir
    }

    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f {
        float x, y, z;
        Vec3f() : x(0), y(0), z(0) {}
        Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}


        Vec3f operator+(const Vec3f &v) const {
            return Vec3f{x + v.x, y + v.y, z + v.z};
        }

        Vec3f operator+(float f) const {
            return Vec3f{x + f, y + f, z + f};
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

        Vec3f operator/(Vec3f v) const {
            return Vec3f{x / v.x, y / v.y, z / v.z};
        }

        Vec3f crossProduct(const Vec3f &v) const {
            return Vec3f{y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x};
        }

        float distance(const Vec3f &v) const {
            return sqrt(pow(x - v.x, 2) + pow(y - v.y, 2) + pow(z - v.z, 2));
        }

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

        bool operator==(const Vec3f &v) {
            return floatEquals(x, v.x) && floatEquals(y, v.y) && floatEquals(z, v.z);
        }

        bool operator>(const Vec3f &v) {
            return x > v.x && y > v.y && z > v.z;
        }


    };

    struct Vec3i {
        int x, y, z;

        Vec3i() {
            x = 0;
            y = 0;
            z = 0;
        }

        Vec3i(int x, int y, int z) : x(x), y(y), z(z) {}

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

    struct Sphere;

    struct Sphere {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct MySphere {
        Vec3f c;
        float r;

        MySphere(Vec3f c, float r) : c(c), r(r) {}
    };

    struct Box {
        Vec3f min, max;

        bool operator==(const Box &b) {
            return min == b.min && max == b.max;
        }

        Vec3f dimensions() {
            return max - min;
        }

        bool is3d() {
            auto dims = dimensions();
            return !floatEquals(dims.x, 0) && !floatEquals(dims.y, 0) && !floatEquals(dims.z, 0);
        }

        bool is2d() {
            auto dims = dimensions();
            if (is3d()) {
                return false;
            }
            if (floatEquals(dims.x, 0)) {
                return !floatEquals(dims.y, 0) && !floatEquals(dims.z, 0);
            }
            if (floatEquals(dims.y, 0)) {
                return !floatEquals(dims.x, 0) && !floatEquals(dims.z, 0);
            }
            if (floatEquals(dims.z, 0)) {
                return !floatEquals(dims.x, 0) && !floatEquals(dims.y, 0);
            }
            return false;
        }

        bool contains(Vec3f point) {
            return point.x >= min.x && point.x <= max.x &&
                   point.y >= min.y && point.y <= max.y &&
                   point.z >= min.z && point.z <= max.z;
        };
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

    struct TriangleVertex {
        int v_id;
        Face *face;
        int material_id;
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

        MySphere getBoundingSphere(Mesh &mesh) {
            Vec3f min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                         std::numeric_limits<float>::max()};
            Vec3f max = {-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
                         -std::numeric_limits<float>::max()};
            for (auto &face: mesh.faces) {
                for (auto &vertexId: {face.v0_id, face.v1_id, face.v2_id}) {
                    auto &vertex = vertex_data[vertexId - 1];
                    for (int i = 0; i < 3; i++) {
                        if (vertex[i] < min[i]) {
                            min[i] = vertex[i];
                        }
                        if (vertex[i] > max[i]) {
                            max[i] = vertex[i];
                        }
                    }
                }
            }
            auto center = (min + max) * 0.5;
            auto radius = center.distance(max);
            return {center, radius};
        }

        Box getBoundingBox(std::vector<int> &v_ids) {
            Vec3f min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
                         std::numeric_limits<float>::max()};
            Vec3f max = {-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
                         -std::numeric_limits<float>::max()};
            for (auto i: v_ids) {
                auto &vertex = vertex_data[i - 1];
                for (int j = 0; j < 3; j++) {
                    if (vertex[j] < min[j]) {
                        min[j] = vertex[j];
                    }
                    if (vertex[j] > max[j]) {
                        max[j] = vertex[j];
                    }
                }
            }
            return {min, max};
        }
    };
}

#endif
