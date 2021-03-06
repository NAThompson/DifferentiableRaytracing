#ifndef DRT_HALLEY_HPP
#define DRT_HALLEY_HPP
#include <utility>
#include <tuple>
#include <cmath>
#include <iostream>
#include <drt/roots.hpp>
#include <drt/vec.hpp>
#include <drt/matrix.hpp>
#include <drt/tensor.hpp>
#include <drt/bounds.hpp>

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

// See Cuyt, "Computational Implementation of the Multivariate Halley Method for Solving Nonlinear Systems of Equations."
template<typename Real, int64_t dimension>
vec<Real, dimension> halley(std::function<std::tuple<vec<Real, dimension>, matrix<Real,dimension,dimension>, tensor<Real,dimension,dimension,dimension>>(vec<Real, dimension>)> f, bounds<Real, dimension> bound, vec<Real, dimension> guess) {
    using std::abs;
    using std::sqrt;
    Real eps = std::numeric_limits<Real>::epsilon();
    auto [v, J, H] = f(guess);
    Real g = squared_norm(v)/2;
#ifdef DEBUG
    std::cerr << std::fixed;
    std::cerr << dimension << " dimensional Halley: guess,        f(guess),           ‖f‖²/2,          δw,        iteration\n";
#endif

    int i = 1;
    do {
        auto a = J.solve(-v);
#ifdef DEBUG
        std::cerr << "\t" << guess << ", " << v << ", " << g << ", " << a << ", " << i << "\n";
#endif
        auto Haa = H(a,a);
        auto b = J.solve(Haa);
        vec<Real,dimension> newguess;
        for (size_t i = 0; i < dimension; ++i) {
            Real denom = a[i] + b[i]/2;
            if (denom == 0) {
                newguess[i] = guess[i] + a[i];
            } else {
                newguess[i] = guess[i] + a[i]*a[i]/(a[i] + b[i]/2);
            }
        }
        if (bound.contains(newguess))
        {
            guess = newguess;
        }
        else
        {
            //std::cerr << __FILE__ << ":" << __LINE__ << " Out of bounds!\n";
            guess = bound.random();
        }

        std::tie(v, J, H) = f(guess);
        g = squared_norm(v)/2;

    // TODO: Some actual analysis on this termination condition:
    } while (g > 100*eps*(abs(v[0]) + abs(v[1])) && i++ < 30);

#ifdef DEBUG
        std::cerr << "\t" << guess << ", " << v << ", " << g << "\n";
#endif
    return guess;
}


}
#endif