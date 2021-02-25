#ifndef DRT_HALLEY_HPP
#define DRT_HALLEY_HPP
#include <utility>
#include <tuple>
#include <cmath>
#include <iostream>
#include <drt/roots.hpp>
#include <drt/vec.hpp>
#include <drt/mat.hpp>

namespace drt {

template<typename Real>
Real bisect_backup(std::function<std::tuple<Real,Real,Real>(Real)> f, Real tmin, Real tmax) {
    auto g = [&](Real x)->Real {
        return std::get<0>(f(x));
    };
    auto [t_min, t_max] = bisect<Real>(g, tmin, tmax);
    return (t_min+t_max)/2;
}



// Newtons method, but more useful for compute graphics since it makes some attempt to recover the minimal t solution.
template<typename Real>
Real halley(std::function<std::tuple<Real,Real,Real>(Real)> f, Real tmin, Real tmax) {
    using std::abs;
    using std::sqrt;
    Real t = tmin;
    Real eps = std::numeric_limits<Real>::epsilon();
    auto [y, dydt, d2ydt2] = f(t);
#ifdef DEBUG
//std::cerr << std::setprecision(std::numeric_limits<Real>::digits10 + 1);
//std::cerr << "Halley iterate:\n";
#endif
    int iterations = 0;
    do {
        if (y == 0) {
            return t;
        }
        if (dydt == 0) {
            // This is a special case of Halley's irrational method:
            // 0 = f(x+h) = f(x) + hf'(x) + h²f''(x)/2 => h² = -2f(x)/f''(x)
            if (d2ydt2 == 0) {
                t = bisect_backup(f, tmin, tmax);
            }
            Real hsq = -2*y/d2ydt2;
            if (hsq < 0) {
                t = bisect_backup(f, tmin, tmax);
            }
            Real h = sqrt(hsq);
            // Always go for the minimal t root:
            t -= h;
            if (t > tmax || t < tmin) {
                t += 2*h;
            }
        }
        else {
            Real h = -2*y*dydt/(2*dydt*dydt - y*d2ydt2);
            t += h;
            if (t < tmin || t > tmax) {
                t = bisect_backup(f, tmin, tmax);
            }
        }
#ifdef DEBUG
        //std::cerr << "\t(t, f(t), f'(t), f''(t), i) = ("  << t << ", " << y << ", "  << dydt << ", " << d2ydt2 << ", " << iterations << ")\n";
#endif
        std::tie(y, dydt, d2ydt2) = f(t);
    } while (abs(y) > eps*abs(t*dydt) && iterations++ < 32);

    std::tie(y, dydt, d2ydt2) = f(t);
    if (abs(y) > eps*abs(t*dydt)) {
        return bisect_backup(f, tmin, tmax);
    }
    return t;
}


}
#endif