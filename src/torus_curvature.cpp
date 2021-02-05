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

using std::make_shared;
using namespace drt;

template<typename Real>
drt::vec<Real, 3> ray_color(const drt::ray<Real>& r, const drt::vec<Real>& background_color, const drt::hittable<Real> & world, int depth) {
    if (depth <= 0) {
        return drt::vec<Real>(0,0,0);
    }
    drt::hit_record<Real> rec;
    if (!world.hit(r, 1000*std::numeric_limits<Real>::epsilon(), std::numeric_limits<Real>::infinity(), rec)) {
        return background_color;
    }
    drt::ray<Real> scattered;
    drt::vec<Real> attenuation;
    drt::vec<Real> emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
        //std::cout << "Returning the emitted because it doesn't scatter; emitted is " << emitted << "\n";
        return emitted;
    }
    drt::vec<Real> color = ray_color(scattered, background_color, world, depth-1);
    // Peter Shirley has an overload for componentwise vector multiplication.
    // I'm avoiding that for now.
    // It does make sense for this application: Independently attenuating each frequency.
    //std::cout << "Emitted color is " << emitted << "\n";
    //std::cout << "Pre-attenuation color is " << color << "\n";
    //std::cout << "Attenuation is " << attenuation << "\n";
    for (size_t i = 0; i < 3; ++i) {
        color[i] *= attenuation[i];
    }
    return emitted + color;
}

int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int>(image_width);
    const int64_t samples_per_pixel = 8;

    auto mat = make_shared<lambertian<Real>>(vec<Real>(1,0,0));
    vec<Real> center(0,0,0);
    auto world = torus<Real>(center, 3.0, 1.0, mat);

    drt::vec<Real> lookfrom(0, 0, -10);
    drt::vec<Real> lookat = center;
    drt::vec<Real> vup(0,1,0);

    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937_64 gen;
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    drt::vec<Real> background(0.01, 0.01, 0.01);
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
