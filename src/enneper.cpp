#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>
#include <drt/lambertian.hpp>
#include <drt/texture.hpp>
#include <drt/diffuse_light.hpp>
#include <drt/aarect.hpp>
#include <drt/color_maps.hpp>
#include <drt/ray_color.hpp>
#include <drt/render_scene.hpp>
#include <drt/enneper.hpp>
#include <drt/dielectric.hpp>
#include <drt/cylinder.hpp>

using std::make_shared;
using std::log;
using namespace drt;

template<typename Real>
hittable_list<Real> enneper_scene() {
    hittable_list<Real> objects;
    Real radius = 4;
    auto enne = make_shared<enneper<Real>>(radius);
    auto mat = make_shared<diffuse_light<Real>>(vec<Real,3>(0.02, 0.4, 0.9));
    objects.add(enne, mat);
    return objects;
}


int main() {
    using Real = double;

    const Real aspect_ratio = 1.0;
    const int64_t image_width = 800;
    const int64_t image_height = static_cast<int64_t>(image_width/aspect_ratio);
    const int64_t samples_per_pixel = 1;

    auto world = enneper_scene<Real>();
    drt::vec<Real> lookat(0,0,0);
    drt::vec<Real> lookfrom(0, -10, 0);
    drt::vec<Real> vup(0,0,1);
    drt::camera<Real> cam(lookfrom, lookat, vup, Real(40), aspect_ratio);
    drt::vec<Real> background(0.0, 0.0, 0.0);
    drt::render_scene<Real>("enneper.png", image_width, image_height, background, cam, world, samples_per_pixel);
}
