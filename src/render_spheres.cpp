#include <iostream>
#include <random>
#include <drt/vec.hpp>
#include <drt/png.hpp>
#include <drt/ray.hpp>
#include <drt/sphere.hpp>
#include <drt/hittable.hpp>
#include <drt/hittable_list.hpp>
#include <drt/camera.hpp>


template<typename Real>
drt::vec<Real, 3> ray_color(const drt::ray<Real>& r, const drt::hittable<Real> & world) {
    drt::hit_record<Real> rec;
    if (world.hit(r, 0, std::numeric_limits<Real>::infinity(), rec)) {
        return 0.5 * (rec.normal + drt::vec<Real, 3>(1,1,1));
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
    world.add(std::make_shared<drt::sphere<Real>>(drt::vec<Real,3>(0,0,-1), 0.5));
    world.add(std::make_shared<drt::sphere<Real>>(drt::vec<Real,3>(0,-100.5,-1), 100));
    drt::camera<Real> cam;
    std::uniform_real_distribution<Real> dis(0,1);
    std::mt19937 gen;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    for (int64_t j = 0; j < image_height; ++j) {
        for (int64_t i = 0; i < image_width; ++i) {
            drt::vec<Real, 3> color(0,0,0);
            for (int64_t s = 0; s < samples_per_pixel; ++s) {
                Real u = (Real(i) + dis(gen)) / (image_width -1);
                Real v = (Real(j) + dis(gen)) / (image_height -1);
                auto r = cam.get_ray(u, v);
                color += ray_color(r, world);
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
