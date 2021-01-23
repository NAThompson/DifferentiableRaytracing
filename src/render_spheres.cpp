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


template<typename Real>
drt::vec<Real, 3> ray_color(const drt::ray<Real>& r, const drt::hittable<Real> & world, int depth) {
    if (depth <= 0) {
        return drt::vec<Real>(0,0,0);
    }
    drt::hit_record<Real> rec;
    if (world.hit(r, std::sqrt(std::numeric_limits<Real>::epsilon()), std::numeric_limits<Real>::infinity(), rec)) {
        drt::ray<Real> scattered;
        drt::vec<Real> attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            drt::vec<Real> color = ray_color(scattered, world, depth-1);
            // Peter Shirley has an overload for componentwise vector multiplication.
            // I'm avoiding that for now.
            // It does make sense for this application: Independently attenuating each frequency.
            for (size_t i = 0; i < 3; ++i) {
                color[i] *= attenuation[i];
            }
            return color;
        }

        return drt::vec<Real>(0,0,0);
    }

    auto n = r.direction();
    drt::normalize(n);
    auto t = (n[1] + 1)/2;
    return (1-t)*drt::vec<Real, 3>(1.0, 1.0, 1.0) + t*drt::vec<Real, 3>(0.5, 0.7, 1.0);
}

int main() {
    using Real = double;

    const auto aspect_ratio = 16.0 / 9.0;
    const int64_t image_width = 1880;
    const int64_t image_height = static_cast<int>(image_width / aspect_ratio);
    const int64_t samples_per_pixel = 16;

    drt::hittable_list<Real> world;

    auto material_ground = std::make_shared<drt::lambertian<Real>>(drt::vec<Real>(0.8, 0.8, 0.0));
    auto material_center = std::make_shared<drt::lambertian<Real>>(drt::vec<Real>(0.1, 0.2, 0.5));
    auto material_left   = std::make_shared<drt::dielectric<Real>>(1.5);
    auto material_right  = std::make_shared<drt::metal<Real>>(drt::vec<Real>(0.8, 0.6, 0.2));

    world.add(std::make_shared<drt::sphere<Real>>(drt::vec<Real>( 0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(std::make_shared<drt::sphere<Real>>(drt::vec<Real>( 0.0,    0.0, -1.0),   0.5, material_center));
    world.add(std::make_shared<drt::sphere<Real>>(drt::vec<Real>(-1.0,    0.0, -1.0),   0.5, material_left));
    world.add(std::make_shared<drt::sphere<Real>>(drt::vec<Real>( 1.0,    0.0, -1.0),   0.5, material_right));

    drt::camera<Real> cam;
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937 gen;
    int max_depth = 5;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    for (int64_t j = 0; j < image_height; ++j) {
        for (int64_t i = 0; i < image_width; ++i) {
            drt::vec<Real, 3> color(0,0,0);
            for (int64_t s = 0; s < samples_per_pixel; ++s) {
                Real u = (Real(i) + dis(gen)) / (image_width -1);
                Real v = (Real(j) + dis(gen)) / (image_height -1);
                auto r = cam.get_ray(u, v);
                color += ray_color(r, world, max_depth);
            }
            auto c = drt::to_8bit_rgba(color/Real(samples_per_pixel));
            int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
        }
    }

    drt::write_png("first.png", img, image_width, image_height);
}
