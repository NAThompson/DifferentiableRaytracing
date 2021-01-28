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

template<typename Real>
drt::hittable_list<Real> simple_light() {
    using std::make_shared;
    using drt::sphere;
    using drt::vec;
    using drt::lambertian;
    using drt::diffuse_light;
    using drt::constant_medium;
    drt::hittable_list<Real> objects;

    auto chtext = make_shared<drt::checker_texture<Real>>(vec<Real>(0.8, 0.0, 1.0), vec<Real>(0.0, 0.8, 1.0));
    objects.add(make_shared<sphere<Real>>(vec<Real>(0,-1000,0), Real(1000), make_shared<drt::lambertian<Real>>(chtext)));

    auto stext = make_shared<drt::solid_color<Real>>(vec<Real>(0.2, 0.4, 0.9));
    objects.add(make_shared<sphere<Real>>(vec<Real>(0,2,0), Real(2), make_shared<drt::metal<Real>>(stext)));

    auto difflight = make_shared<diffuse_light<Real>>(vec<Real>(4,4,4));
    objects.add(make_shared<drt::xy_rect<Real>>(Real(3), Real(5), Real(1), Real(3), Real(-2), difflight));

    return objects;
}

template<typename Real>
drt::vec<Real, 3> ray_color(const drt::ray<Real>& r, const drt::vec<Real>& background_color, const drt::hittable<Real> & world, int depth) {
    if (depth <= 0) {
        return drt::vec<Real>(0,0,0);
    }
    drt::hit_record<Real> rec;
    if (!world.hit(r, std::sqrt(std::numeric_limits<Real>::epsilon()), std::numeric_limits<Real>::infinity(), rec)) {
        return background_color;
    }
    drt::ray<Real> scattered;
    drt::vec<Real> attenuation;
    drt::vec<Real> emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
        return emitted;
    }
    drt::vec<Real> color = ray_color(scattered, background_color, world, depth-1);
    // Peter Shirley has an overload for componentwise vector multiplication.
    // I'm avoiding that for now.
    // It does make sense for this application: Independently attenuating each frequency.
    for (size_t i = 0; i < 3; ++i) {
        color[i] *= attenuation[i];
    }
    return emitted + color;
}

int main() {
    using Real = float;

    const Real aspect_ratio = 1.6;
    const int64_t image_width = 1800;
    const int64_t image_height = static_cast<int>(image_width / aspect_ratio);
    const int64_t samples_per_pixel = 16;

    drt::hittable_list<Real> world1 = simple_light<Real>();
    drt::bvh_node<Real> world(world1);

    drt::vec<Real> lookfrom(26,3,6);
    drt::vec<Real> lookat(0,2,0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(20), aspect_ratio);
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937 gen;
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    drt::vec<Real> background(0.01, 0.02, 0.01);
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

    drt::write_png("lights.png", img, image_width, image_height);
}