#include <complex>
#include <tuple>
#include <iostream>
#include <chrono>
#include <drt/png.hpp>
#include <drt/color_maps.hpp>
#include <drt/progress_bar.hpp>


template<typename Real>
auto fifth_roots(std::complex<Real> z) {
    std::complex<Real> v = std::pow(z,4);
    std::complex<Real> dw = Real(5)*v;
    std::complex<Real> w = v*z - Real(1);
    return std::make_pair(w, dw);
}

template<typename Real>
auto g(std::complex<Real> z) {
    std::complex<Real> z2 = z*z;
    std::complex<Real> z3 = z*z2;
    std::complex<Real> z4 = z2*z2;
    std::complex<Real> w = z4*(z4 + Real(15)) - Real(16);
    std::complex<Real> dw = Real(4)*z3*(Real(2)*z4 + Real(15));
    return std::make_pair(w, dw);
}

template<typename Real>
auto fifth_roots_halley(std::complex<Real> z) {
    std::complex<Real> zcb = std::pow(z,3);
    std::complex<Real> dw = Real(5)*zcb*z;
    std::complex<Real> w = zcb*z*z - Real(1);
    std::complex<Real> ddw = Real(20)*zcb;
    return std::make_tuple(w, dw, ddw);
}

template<typename Real>
std::complex<Real> complex_newton(std::function<std::pair<std::complex<Real>,std::complex<Real>>(std::complex<Real>)> f, std::complex<Real> z) {
    // f(x(1+e)) = f(x) + exf'(x)
    bool close = false;
    do {
        auto [y, dy] = f(z);
        z -= y/dy;
        close = (abs(y) <= 1.4*std::numeric_limits<Real>::epsilon()*abs(z*dy));
    } while(!close);
    return z;
}

template<typename Real>
std::complex<Real> complex_halley(std::function<std::tuple<std::complex<Real>,std::complex<Real>,std::complex<Real>>(std::complex<Real>)> f, std::complex<Real> z) {
    // f(x(1+e)) = f(x) + exf'(x)
    bool close = false;
    do {
        auto [y, dy, ddy] = f(z);
        z -= (Real(2)*y*dy/(Real(2)*dy*dy - y*ddy));
        close = (abs(y) <= std::numeric_limits<Real>::epsilon()*abs(z*dy)/2);
    } while(!close);
    return z;
}

template<typename Real>
class plane_pixel_map
{
public:
    plane_pixel_map(int64_t image_width, int64_t image_height, Real xmin, Real ymin)
    {
        image_width_ = image_width;
        image_height_ = image_height;
        xmin_ = xmin;
        ymin_ = ymin;
    }

    std::complex<Real> to_complex(int64_t i, int64_t j) const {
        Real x = xmin_ + 2*abs(xmin_)*Real(i)/Real(image_width_ - 1);
        Real y = ymin_ + 2*abs(ymin_)*Real(j)/Real(image_height_ - 1);
        return std::complex<Real>(x,y);
    }

    std::pair<int64_t, int64_t> to_pixel(std::complex<Real> z) const {
        Real x = z.real();
        Real y = z.imag();
        Real ii = (image_width_ - 1)*(x - xmin_)/(2*abs(xmin_));
        Real jj = (image_height_ - 1)*(y - ymin_)/(2*abs(ymin_));

        return std::make_pair(std::round(ii), std::round(jj));
    }

private:
    int64_t image_width_;
    int64_t image_height_;
    Real xmin_;
    Real ymin_;
};

int main() {
    using Real = long double;
    int64_t image_width = 4096;
    int64_t image_height = 4096;
    std::vector<uint8_t> img(4*image_width*image_height, 0);
    plane_pixel_map<Real> map(image_width, image_height, Real(-2), Real(-2));
    for (int64_t j = 0; j < image_height; ++j) {
        drt::display_progress(Real(j)/Real(image_height));
        for (int64_t i = 0; i < image_width; ++i) {
            std::complex<Real> z0 = map.to_complex(i,j);
            //auto rt = complex_halley<Real>(fifth_roots_halley<Real>, z0);
            auto rt = complex_newton<Real>(g<Real>, z0);
            // The root is one of exp(2πij/5). Therefore is can be classified by angle.
            Real theta = atan2(rt.imag(), rt.real());
            // Now theta in [-π,π]. Get it into [0,2π]:
            if (theta < 0) {
                theta += 2*M_PI;
            }
            theta /= 2*M_PI;
            if (std::isnan(theta)) {
                std::cerr << "Theta is a nan!\n";
            }
            auto c = drt::to_8bit_rgba(drt::smooth_cool_warm(theta));
            int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
            img[idx + 0] = c[0];
            img[idx + 1] = c[1];
            img[idx + 2] = c[2];
            img[idx + 3] = c[3];
        }
    }

    std::array<std::complex<Real>, 8> roots;
    roots[0] = -Real(1);
    roots[1] = Real(1);
    roots[2] = {Real(0), Real(1)};
    roots[3] = {Real(0), -Real(1)};
    roots[4] = {sqrt(Real(2)), sqrt(Real(2))};
    roots[5] = {sqrt(Real(2)), -sqrt(Real(2))};
    roots[6] = {-sqrt(Real(2)), -sqrt(Real(2))};
    roots[7] = {-sqrt(Real(2)), sqrt(Real(2))};
    for (int64_t k = 0; k < 8; ++k) {
        auto [ic, jc] = map.to_pixel(roots[k]);

        int64_t r = 7;
        for (int64_t i = ic - r; i < ic + r; ++i) {
            for (int64_t j = jc - r; j < jc + r; ++j) {
                if ((i-ic)*(i-ic) + (j-jc)*(j-jc) > r*r) {
                    continue;
                }
                int64_t idx = 4 * image_width * (image_height - 1 - j) + 4 * i;
                img[idx + 0] = 0;
                img[idx + 1] = 0;
                img[idx + 2] = 0;
                img[idx + 3] = 0xff;
            }
        }
    }

    drt::write_png("newton_fractal.png", img, image_width, image_height);
}
