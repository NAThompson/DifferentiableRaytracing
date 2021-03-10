#ifndef DRT_NEWTON_HPP
#define DRT_NEWTON_HPP
#include <utility>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <drt/roots.hpp>
#include <drt/vec.hpp>
#include <drt/matrix.hpp>
#include <drt/bounds.hpp>

#ifndef DRT_NEWTON_DEBUG
//#define DRT_NEWTON_DEBUG
#endif

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
#ifdef NEWTON_DEBUG
//std::cerr << std::setprecision(std::numeric_limits<Real>::digits10 + 1);
//std::cerr << "Newton iterate:\n";
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
#ifdef NEWTON_DEBUG
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
    return t;
}


// See Numerical Recipes, section 9.7: Globally Convergent Methods for Nonlinear Systems of Equations.
template<typename Real, int64_t dimension>
vec<Real, dimension> newton(std::function<std::pair<vec<Real, dimension>, matrix<Real,dimension,dimension>>(vec<Real, dimension>)> f,
                            bounds<Real,dimension> const & bound, vec<Real,dimension> guess, int64_t max_iterations = 16) {
    using std::abs;
    using std::sqrt;
    Real eps = std::numeric_limits<Real>::epsilon();
#ifdef DRT_NEWTON_DEBUG
    //std::cerr << std::fixed;
    std::cerr << std::setprecision(std::numeric_limits<Real>::digits10+1);
    std::cerr << "2D Newton: x,        f(x),              ‖f‖^2/2,          δw,        iteration\n";
#endif
    auto [v, J] = f(guess);
    Real g = squared_norm(v)/2;
    int i = 1;
    Real expected_residual = 0;
    int left_domain = 0;
    do {
        auto dw = J.solve(-v);
#ifdef DRT_NEWTON_DEBUG
        std::cerr << "\t" << guess << ", " << v << ", " << g << ", " << dw << ", " << i << "\n";
#endif
        vec<Real, dimension> newguess = guess + dw;

        if (!bound.contains(newguess)) {
            // This assert fails for trivial reasons: A single ulp oob.
            //assert(bound.contains(guess));
            Real lambda_max = 1;
            for (int64_t i = 0; i < dimension; ++i) {
                if (newguess[i] < bound[i].first) {
                    lambda_max = std::min((bound[i].first - guess[i])/dw[i], lambda_max);
                }
                if (newguess[i] > bound[i].second) {
                    lambda_max = std::min((bound[i].second - guess[i])/dw[i], lambda_max);
                }
            }

            // I'm not sure about this . . .
            // For sure, if lambda_max <= 0, this has to be done.
            // This is just a fudge factor saying "we're really close to the boundary, and Newton's method wants us to move out."
            // "Find a new starting location."
            if (lambda_max <= 0.1) {
                // Oh boy how bad is this.
                guess = bound.random();
                assert(bound.contains(guess));
                #ifdef DRT_NEWTON_DEBUG
                std::cerr << "\tRandomizing as update is out of domain; random choice: x = " << guess << ".\n";
                #endif
                ++left_domain;
                if (left_domain > 6) {
                    return vec<Real,dimension>(special_vec::NaNs);
                }
                std::tie(v, J) = f(guess);
                g = squared_norm(v)/2;
                continue;
            }

            assert(lambda_max < 1);
            #ifdef DRT_NEWTON_DEBUG
            std::cerr << "\tNewton step leaves domain; taking step λ = " << lambda_max << "\n";
            #endif
            dw *= lambda_max;
            newguess = guess + dw;
        }

        auto [vnew, Jnew] = f(newguess);
        Real gnew = squared_norm(vnew)/2;
        if (gnew > g) {
            // Backtrack: Numerical Recipes 9.7.11:
            Real g1 = gnew;
            Real g0 = g;
            Real gprime0 = dot(J*dw, v);
            Real lambda = -gprime0/(2*(g1 - g0 - gprime0));
            #ifdef DRT_NEWTON_DEBUG
            std::cerr << "\tBacktracking. Chosen λ = " << lambda << "\n";
            #endif
            newguess = guess + lambda*dw;
            std::tie(vnew, Jnew) = f(newguess);
            gnew = squared_norm(vnew)/2;
            if(gnew > g) {
                #ifdef DRT_NEWTON_DEBUG
                std::cerr << "\tBacktracking made things worse. Randomizing.\n";
                #endif
                newguess = bound.random();
                assert(bound.contains(newguess));
                std::tie(vnew, Jnew) = f(newguess);
                gnew = squared_norm(v)/2;
            }
        }

        g = gnew;
        guess = newguess;
        v = vnew;
        J = Jnew;
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
        #ifdef DRT_NEWTON_DEBUG
        std::cerr << "\tResidual = " << max_norm(v) << ", expected residual = " << expected_residual << "\n";
        #endif
    } while (max_norm(v) > expected_residual && i++ < max_iterations);

#ifdef DRT_NEWTON_DEBUG
        std::cerr << "\t" << guess << ", " << v << ", " << g << "\n";
#endif
    if (max_norm(v) <= expected_residual) {
        return guess;
    }
    else {
        return vec<Real,dimension>(special_vec::NaNs);
    }
    return guess;
}

}
#endif