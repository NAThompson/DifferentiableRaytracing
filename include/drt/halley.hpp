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

#ifndef DRT_HALLEY_DEBUG
//#define DRT_HALLEY_DEBUG
#endif

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
#ifdef DRT_HALLEYDEBUG
std::cerr << std::setprecision(std::numeric_limits<Real>::digits10 + 1);
std::cerr << "Halley iterate:\n";
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
#ifdef DRT_HALLEY_DEBUG
        std::cerr << "\t(t, f(t), f'(t), f''(t), i) = ("  << t << ", " << y << ", "  << dydt << ", " << d2ydt2 << ", " << iterations << ")\n";
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
// TODO: Also implement backtracking.
template<typename Real, int64_t dimension>
vec<Real, dimension> halley(std::function<std::tuple<vec<Real, dimension>, matrix<Real,dimension,dimension>, tensor<Real,dimension,dimension,dimension>>(vec<Real, dimension>)> f,
                            bounds<Real, dimension> bound,
                            vec<Real, dimension> guess,
                            int max_iterations=16) {
    using std::abs;
    using std::sqrt;
    Real eps = std::numeric_limits<Real>::epsilon();
    auto [v, J, H] = f(guess);
    Real g = squared_norm(v)/2;
#ifdef DRT_HALLEY_DEBUG
    std::cerr << std::fixed;
    std::cerr << dimension << " dimensional Halley: guess,        f(guess),           ‖f‖²/2,          δw,        iteration\n";
#endif

    int i = 1;
    int randomizations = 0;
    Real expected_residual = 0;
    do {
        auto a = J.solve(-v);
#ifdef DRT_HALLEY_DEBUG
        std::cerr << "\t" << guess << ", " << v << ", " << g << ", " << a << ", " << i << "\n";
#endif
        auto Haa = H(a,a);
        auto b = J.solve(Haa);
        vec<Real,dimension> newguess;
        bool randomized = false;
        for (size_t i = 0; i < dimension; ++i) {
            Real denom = a[i] + b[i]/2;
            if (denom == 0) {
                newguess[i] = guess[i] + a[i];
            } else {
                newguess[i] = guess[i] + a[i]*a[i]/(a[i] + b[i]/2);
            }
        }
        if (!bound.contains(newguess))
        {
            randomized = true;
            if (randomizations++ < 5) {
                guess = bound.random();
            }
            else {
                vec<Real, dimension> nans;
                for (int64_t i = 0; i < dimension; ++i) {
                    nans[i] = std::numeric_limits<Real>::quiet_NaN();
                }
                return nans;
            }
        } else {
            guess = newguess;
        }

        std::tie(v, J, H) = f(guess);
        Real gnew = squared_norm(v)/2;
        if (gnew > g && !randomized) {
            #ifdef DRT_HALLEY_DEBUG
            std::cerr << "\tQuadratic form increased even though we did not randomize. Need to implement backtracking.\n";
            #endif
        }
        g = gnew;
        expected_residual = 0;
        for (int64_t i = 0; i < dimension; ++i) {
            Real row_residual = 0;
            for (int64_t j = 0; j < dimension; ++j) {
                row_residual += std::abs(J(i,j)*guess[j]);
            }
            // eps/2 is a little harsh, right?
            row_residual *= eps;
            if (row_residual > expected_residual) {
                expected_residual = row_residual;
            }
        }

    } while (max_norm(v) > expected_residual && i++ < max_iterations);

#ifdef DRT_HALLEY_DEBUG
        std::cerr << "\t" << guess << ", " << v << ", " << g << "\n";
#endif

    if (max_norm(v) <= expected_residual) {
        return guess;
    }
    else {
        vec<Real, dimension> nans;
        for (int64_t i = 0; i < dimension; ++i) {
            nans[i] = std::numeric_limits<Real>::quiet_NaN();
        }
        return nans;
    }
}


}
#endif
