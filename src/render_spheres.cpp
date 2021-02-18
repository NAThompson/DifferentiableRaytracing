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
#include <drt/render_scene.hpp>

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

    auto checker = std::make_shared<drt::checker_texture<Real>>(vec<Real>(0.2, 0.3, 0.1), vec<Real>(0.9, 0.9, 0.9));
    world.add(std::make_shared<drt::sphere<Real>>(vec<Real>(0,-1000,0), 1000), std::make_shared<lambertian<Real>>(checker));

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
                    world.add(make_shared<sphere<Real>>(center, 0.2), sphere_material);
                } else if (choose_mat < 0.75) {
                    // metal
                    vec<Real> albedo;
                    albedo[0] = dis(gen)/2 + 0.5;
                    albedo[1] = dis(gen)/2 + 0.5;
                    albedo[2] = dis(gen)/2 + 0.5;
                    sphere_material = make_shared<metal<Real>>(albedo);
                    world.add(make_shared<sphere<Real>>(center, 0.2), sphere_material);
                } else {
                    // glass
                    sphere_material = make_shared<dielectric<Real>>(1.5);
                    world.add(make_shared<sphere<Real>>(center, 0.2), sphere_material);
                }
            }
        }
    }

    auto material1 = make_shared<dielectric<Real>>(1.5);
    world.add(make_shared<sphere<Real>>(vec<Real>(0, 1, 0), 1.0), material1);

    auto material2 = make_shared<lambertian<Real>>(vec<Real>(0.4, 0.2, 0.1));
    world.add(make_shared<sphere<Real>>(vec<Real>(-4, 1, 0), 1.0), material2);

    auto material3 = make_shared<metal<Real>>(vec<Real>(0.7, 0.6, 0.5));
    world.add(make_shared<sphere<Real>>(vec<Real>(4, 1, 0), 1.0), material3);

    return world;
}

int main() {
    using Real = float;

    const Real aspect_ratio = 1.6;
    const int64_t image_width = 1800;
    const int64_t image_height = static_cast<int>(image_width / aspect_ratio);
    const int64_t samples_per_pixel = 128;

    drt::hittable_list<Real> world = random_scene<Real>();

    drt::vec<Real> lookfrom(13,2,3);
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(20), aspect_ratio);
    drt::vec<Real> background(0.70, 0.80, 1.00);
    drt::render_scene<Real>("spheres.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
