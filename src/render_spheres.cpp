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
drt::hittable_list<Real> random_scene() {
    std::uniform_real_distribution<Real> dis(0,1);
    std::random_device gen;
    using std::make_shared;
    using drt::lambertian;
    using drt::vec;
    using drt::sphere;
    using drt::material;
    using drt::metal;
    using drt::dielectric;

    drt::hittable_list<Real> world;

    auto ground_material = make_shared<lambertian<Real>>(vec<Real>(0.5, 0.5, 0.5));
    world.add(make_shared<sphere<Real>>(vec<Real>(0,-1000,0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = dis(gen);
            vec<Real> center(a + 0.9*dis(gen), 0.2, b + 0.9*dis(gen));

            if (norm(center - vec<Real>(4, 0.2, 0)) > 0.9) {
                std::shared_ptr<material<Real>> sphere_material;

                if (choose_mat < 0.5) {
                    // diffuse
                    vec<Real> albedo;
                    albedo[0] = dis(gen)*dis(gen);
                    albedo[1] = dis(gen)*dis(gen);
                    albedo[2] = dis(gen)*dis(gen);
                    sphere_material = make_shared<lambertian<Real>>(albedo);
                    world.add(make_shared<sphere<Real>>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.75) {
                    // metal
                    vec<Real> albedo;
                    albedo[0] = dis(gen)/2 + 0.5;
                    albedo[1] = dis(gen)/2 + 0.5;
                    albedo[2] = dis(gen)/2 + 0.5;
                    sphere_material = make_shared<metal<Real>>(albedo);
                    world.add(make_shared<sphere<Real>>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric<Real>>(1.5);
                    world.add(make_shared<sphere<Real>>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric<Real>>(1.5);
    world.add(make_shared<sphere<Real>>(vec<Real>(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian<Real>>(vec<Real>(0.4, 0.2, 0.1));
    world.add(make_shared<sphere<Real>>(vec<Real>(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal<Real>>(vec<Real>(0.7, 0.6, 0.5));
    world.add(make_shared<sphere<Real>>(vec<Real>(4, 1, 0), 1.0, material3));

    return world;
}

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
    const int64_t samples_per_pixel = 32;

    drt::hittable_list<Real> world = random_scene<Real>();

    drt::vec<Real> lookfrom(13,2,3);
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, 20.0, aspect_ratio);
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937 gen;
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    for (int64_t j = 0; j < image_height; ++j) {
        std::cerr << j << "/" << image_height << "\r";
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
