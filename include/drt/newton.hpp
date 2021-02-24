#ifndef DRT_NEWTON_HPP
#define DRT_NEWTON_HPP
#include <utility>
#include <cmath>
#include <iostream>
#include <drt/roots.hpp>
#include <drt/vec.hpp>
#include <drt/mat.hpp>

namespace drt {

template<typename Real>
Real bisect_backup(std::function<std::pair<Real,Real>(Real)> f, Real tmin, Real tmax) {
    auto g = [&](Real x)->Real {
        return f(x).first;
    };
    auto [t_min, t_max] = bisect<Real>(g, tmin, tmax);
    return (t_min+t_max)/2;
}

// Newtons method, but more useful for compute graphics since it makes some attempt to recover the minimal t solution.
template<typename Real>
Real newton(std::function<std::pair<Real,Real>(Real)> f, Real tmin, Real tmax) {
    using std::abs;
    using std::sqrt;
    Real t = tmin;
    Real eps = std::numeric_limits<Real>::epsilon();
    auto [y, dydt] = f(t);
#ifdef DEBUG
std::cerr << std::setprecision(std::numeric_limits<Real>::digits10 + 1);
std::cerr << "Newton iterate:\n";
#endif
    int iterations = 0;
    do {
        if (y == 0) {
            return t;
        }
        // I used Mueller's method in boost...that's probably better . . .
        if (dydt == 0) {
            t = bisect_backup(f, tmin, tmax);
        }
        else {
            t -= y/dydt;
            if (t < tmin || t > tmax) {
                t = bisect_backup(f, tmin, tmax);
            }
        }
#ifdef DEBUG
        std::cerr << "\t(t, f(t), f'(t), i) = ("  << t << ", " << y << ", "  << dydt << ", " << iterations << ")\n";
#endif
        std::tie(y, dydt) = f(t);
    // This termination is a bit weird. t and abs(y) are not sync'd, intentionally.
    // If it's converging, then the residual is ~cond*eps.
    // If it's not converging well, there's still a good chance it'll terminate.
    } while (abs(y) > eps*abs(t*dydt) && iterations++ < 32);

    std::tie(y, dydt) = f(t);
    if (abs(y) > eps*abs(t*dydt)) {
        return bisect_backup(f, tmin, tmax);
    }
    return t - y/dydt;
}


// See Numerical Recipes, section 9.7: Globally Convergent Methods for Nonlinear Systems of Equations.
template<typename Real>
vec<Real, 2> newton(std::function<std::pair<vec<Real, 2>, mat<Real,2,2>>(Real, Real)> f, Real tmin, Real tmax, Real umin, Real umax) {
    using std::abs;
    using std::sqrt;
    Real t = tmin;
    Real rteps = sqrt(std::numeric_limits<Real>::epsilon());
    auto [v, J] = f(t);
#ifdef DEBUG
    std::cerr << "2D Newton iterate:\n";
#endif

    do {
        Real det = determinant(J);
        if (det == 0) {
#ifdef DEBUG
            std::cerr << "\tNewton's method encountered a zero derivative; moreover bisection failed to find a region of changing sign.\n";
#endif
            return std::numeric_limits<Real>::quiet_NaN();
        }

#ifdef DEBUG
        //std::cerr << "\tResidual = " << y << ", root guess = " << t << ", dydt = " << dydt << "\n";
#endif
    } while (norm(v) > rteps*abs(1));

    return vec<Real,2>(t, umin);
}

}
#endif