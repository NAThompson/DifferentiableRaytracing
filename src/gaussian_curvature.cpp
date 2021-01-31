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

template<typename Real>
drt::hittable_list<Real> ellipsoid_scene() {
    using std::make_shared;
    using drt::lambertian;
    using drt::vec;
    using drt::sphere;
    using drt::dielectric;
    using drt::constant_medium;
    using drt::metal;
    using drt::ellipsoid;

    drt::hittable_list<Real> objects;
    auto light = make_shared<drt::diffuse_light<Real>>(vec<Real>(3, 3, 3));
    objects.add(make_shared<drt::xz_rect<Real>>(100, 423, 100, 412, 700, light));

    Real scale = 3.0;
    Real a = Real(scale*40);
    Real b = Real(scale*70);
    Real c = Real(scale*60);
    vec<Real> center = vec<Real>(260, 250, 45);
    std::function<vec<Real>(Real, Real, const vec<Real> &)> gaussian_curvature = [a,b,c, center]([[maybe_unused]] Real u, [[maybe_unused]] Real v, vec<Real> const & p) -> vec<Real> {
        Real asq = a*a;
        Real bsq = b*b;
        Real csq = c*c;
        // https://mathworld.wolfram.com/Ellipsoid.html, equation 14:
        Real numerator = asq*bsq*bsq*bsq*csq*csq*csq;
        Real sqrt_denom = csq*csq*bsq*bsq + csq*csq*(asq - bsq)*(p[1] -center[1])*(p[1] -center[1]) + bsq*bsq*(asq-csq)*(p[2]-center[2])*(p[2]-center[2]);
        Real kappa = numerator/(sqrt_denom*sqrt_denom);
        Real kappa_min = std::min({csq/(asq*bsq), asq/(csq*bsq), bsq/(asq*csq)});
        Real kappa_max = std::max({csq/(asq*bsq), asq/(csq*bsq), bsq/(asq*csq)});
        if (kappa < kappa_min) {
            std::cout << "kappa is small = " << kappa << ", kappa/kappa_min = " << kappa/kappa_min << "\n";
        }
        if (kappa > kappa_max) {
            std::cout << "kappa is big = " << kappa << "\n";
        }

        Real scalar = (kappa - kappa_min)/(kappa_max - kappa_min);
        vec<Real, 3> w = drt::plasma(scalar);
        return w;
    };

    auto ptr = std::make_shared<decltype(gaussian_curvature)>(gaussian_curvature);
    auto ltext = drt::lambda_texture(ptr);
    auto ltext_ptr = std::make_shared<decltype(ltext)>(ltext);

    auto mat = make_shared<metal<Real>>(ltext_ptr);
    auto boundary = make_shared<ellipsoid<Real>>(center, a, b, c, mat);
    objects.add(boundary);
    return objects;
}

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
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 2880;
    const int64_t image_height = static_cast<int>(image_width / 1.6);
    const int64_t samples_per_pixel = 32;

    drt::hittable_list<Real> world1 = ellipsoid_scene<Real>();
    drt::bvh_node<Real> world(world1);

    drt::vec<Real> lookfrom(359, 278, -600);
    drt::vec<Real> lookat(260, 250,45);
    drt::vec<Real> vup(0,1,0);

    drt::camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937_64 gen;
    int max_depth = 8;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    //drt::vec<Real> background(1.0, 1.0, 1.0);
    drt::vec<Real> background(0.8, 0.8, 0.8);
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

    drt::write_png("gaussian_curvature.png", img, image_width, image_height);
}
