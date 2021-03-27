#include <complex>
#include <tuple>
#include <iostream>
#include <chrono>
#include <drt/png.hpp>
#include <drt/color_maps.hpp>
#include <boost/math/tools/roots.hpp>


template<typename Real>
auto fifth_roots(std::complex<Real> z) {
    std::complex<Real> v = std::pow(z,4);
    std::complex<Real> dw = Real(5)*v;
    std::complex<Real> w = v*z - Real(1);
    return std::make_pair(w, dw);
}


template<typename Real>
std::complex<Real> complex_newton(std::function<std::complex<Real>(std::complex<Real>)> f, std::complex<Real> z) {
    // f(x(1+e)) = f(x) + exf'(x)
    auto [y, dy] = f(z0);
    // .. . .
    return z;
}

template<typename Real>
auto fifth_roots_halley(std::complex<Real> z) {
    std::complex<Real> v = std::pow(z,3);
    std::complex<Real> dw = Real(5)*v;
    std::complex<Real> w = v*z - Real(1);
    return std::make_pair(w, dw);
}

template<typename Real>
std::array<uint8_t, 4> angle_to_color(Real angle) {

    // angle should = 1.25663704*i
    if (angle < 0.1 && angle > -0.1) {
        return {0x51, 0x3b, 0x56, 0xff};
    }
    else if (angle < 1.257 && angle > 1.256) {
        return {0x52, 0x51, 0x74, 0xff};
    }
    else if (angle < 2.52 && angle > 2.51) {
        return {0x34, 0x8a, 0xa7, 0xff};
    }
    else if (angle < 3.77 - 2*M_PI && angle > 3.76 - 2*M_PI) {
        return {0x35, 0x8a, 0xa7, 0xff};
    }
    else if (angle < 5.03 - 2*M_PI && angle > 5.02 - 2*M_PI) {
        return {0x5d, 0xd3, 0x9e, 0xff};
    }
    std::cerr << "Unrecognized angle! angle = " << angle << "\n";
    return {0x00, 0x00, 0x00, 0xff};
}

int main() {
    using Real = long double;
    int64_t image_width = 1024;
    int64_t image_height = 1024;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    for (int64_t j = 0; j < image_height; ++j) {
        for (int64_t i = 0; i < image_width; ++i) {
            Real x = -3 + 6*Real(i)/Real(image_width -1);
            Real y = -3 + 6*Real(j)/Real(image_height -1);
            std::complex<Real> z0(x,y);
            auto rt = boost::math::tools::complex_newton(fifth_roots<Real>, z0, 32);
            // The root is one of exp(2πij/5). Therefore is can be classified by angle.
            Real theta = atan2(rt.imag(), rt.real());
            // Now theta in [-π,π].
            // auto c = angle_to_color(theta);
            if (theta < 0) {
                theta += 2*M_PI;
            }
            theta /= 2*M_PI;
            if (std::isnan(theta)) {
                std::cerr << "Theta is a nan!\n";
            }
            auto c = drt::to_8bit_rgba(drt::viridis(theta));
            int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
        }
    }

    drt::write_png("newton_fractal.png", img, image_width, image_height);

}
