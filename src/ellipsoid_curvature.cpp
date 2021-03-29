#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/png.hpp>
#include <drt/ray.hpp>
#include <drt/aarect.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/material.hpp>
#include <drt/lambertian.hpp>
#include <drt/texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/ellipsoid.hpp>
#include <drt/color_maps.hpp>
#include <drt/render_scene.hpp>
#include <drt/disk.hpp>
#include <drt/dielectric.hpp>
#include <drt/metal.hpp>

#include <omp.h>
using std::make_shared;
using namespace drt;

template<typename Real>
hittable_list<Real> ellipsoid_scene() {
    hittable_list<Real> objects;

    vec<Real> scales(1.95, 1.5, 2.5);
    vec<Real> center = vec<Real>(0, 0, 0);
    auto ell = make_shared<ellipsoid<Real>>(center, scales);
    auto bounds = ell->gaussian_curvature_bounds();
    Real kappa_min = bounds.first;
    Real kappa_max = bounds.second;

    auto gaussian_curvature = [=](hit_record<Real> const & hr) {
        Real kappa = hr.gaussian_curvature();
        Real scalar = (kappa - kappa_min)/(kappa_max - kappa_min);
        return viridis(scalar);
    };

    auto texture = make_shared<lambda_texture<Real>>(gaussian_curvature);
    auto mat = make_shared<lambertian<Real>>(texture);
    //auto mat = make_shared<dielectric<Real>>(1.0, texture);
    //auto mat = make_shared<metal<Real>>(texture);
    objects.add(ell, mat);
    return objects;
}

int main() {
    using Real = double;

    const Real aspect_ratio = 1.6;
    const int64_t image_width = 1200;
    const int64_t image_height = static_cast<int>(image_width / aspect_ratio);
    const int64_t samples_per_pixel = 1024;

    vec<Real> lookat(0, 0, 0);
    vec<Real> vup(0,1,0);
    vec<Real> background(0.0, 0.0, 0.0);
    auto mat = make_shared<diffuse_light<Real>>(vec<Real>(1,1,1));
    #pragma omp parallel for
    for (int64_t i = 0; i < 360; i += 1)
    {
        hittable_list<Real> world = ellipsoid_scene<Real>();
        Real theta = 2*M_PI*Real(i)/Real(360);
        Real r = 7;
        vec<Real> lookfrom(r*cos(theta), 0, r*sin(theta));
        camera cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
        auto [o, d] = cam.backlight();
        auto disk_ptr = make_shared<disk<Real>>(10.0, o, d);
        world.add(disk_ptr, mat);
        std::string filename = "ellipsoid_curvature_";
        if (i < 10) {
            filename += "00";
        }
        else if (i < 100) {
            filename += "0";
        }
        filename += std::to_string(i);
        filename += ".png";
        render_scene<Real>(filename, image_width, image_height, background, cam, world, samples_per_pixel);
    }
}
