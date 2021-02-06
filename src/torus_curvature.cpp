#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/png.hpp>
#include <drt/ray.hpp>
#include <drt/sphere.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/material.hpp>
#include <drt/lambertian.hpp>
#include <drt/dielectric.hpp>
#include <drt/metal.hpp>
#include <drt/aabb.hpp>
#include <drt/bvh.hpp>
#include <drt/texture.hpp>
#include <drt/checker_texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/aarect.hpp>
#include <drt/constant_medium.hpp>
#include <drt/box.hpp>
#include <drt/ellipsoid.hpp>
#include <drt/color_maps.hpp>
#include <drt/torus.hpp>
#include <drt/ray_color.hpp>

using std::make_shared;
using namespace drt;

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int>(image_width);
    const int64_t samples_per_pixel = 8;

    auto mat = make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> center(0,0,0);
    auto world = torus<Real>(center, 3.0, 1.0, mat);

    drt::vec<Real> lookfrom(0, -10, -15);
    drt::vec<Real> lookat = center;
    drt::vec<Real> vup(0,1,0);

    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937_64 gen(12345);
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    drt::vec<Real> background(1.0, 1.0, 1.0);
    for (int64_t j = 0; j < image_height; ++j) {
        std::cerr << j << "/" << image_height << "\r";
        for (int64_t i = 0; i < image_width; ++i) {
            drt::vec<Real, 3> color(0,0,0);
            for (int64_t s = 0; s < samples_per_pixel; ++s) {
                Real u = (Real(i) + dis(gen)) / (image_width -1);
                Real v = (Real(j) + dis(gen)) / (image_height -1);
                auto r = cam.get_ray(u, v);
                color += ray_color(r, background, world, max_depth);
            }
            auto c = drt::to_8bit_rgba(color/Real(samples_per_pixel));
            int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
        }
    }

    drt::write_png("torus_curvature.png", img, image_width, image_height);
}
