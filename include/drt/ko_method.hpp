#ifndef DRT_KO_METHOD_HPP
#define DRT_KO_METHOD_HPP
#include <drt/vec.hpp>
#include <drt/ray.hpp>
#include <drt/hittable.hpp>

namespace drt {

#define DEBUG_KO_METHOD 1

template<typename Real>
vec<Real,3> ko_method(hittable<Real>& h, ray<Real> const & r, Real u0, Real v0, bounds<Real,3> const & bound)
{
    vec<Real> O = r.origin();
    vec<Real> D = r.direction();
    matrix<Real,3,3> J;
    tensor<Real, 3, 2, 2> H;
    vec<Real> P = h(u0,v0, J, H);

    #if DEBUG_KO_METHOD
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    std::cout << "Our initial guess is (u₀,v₀) = (" << u0 << ", " << v0 << ")\n";
    std::cout << "σ(" << u0 << ", " << v0 << ") = " << P << "\n";
    #endif

    vec<Real> PmO = P - O;
    vec<Real> plane_normal = cross(PmO, D);

    // TODO: What is the real condition that these two fall on a line using the rounding model of floating point arithmetic?
    if (squared_norm(plane_normal) < std::numeric_limits<Real>::epsilon()*std::numeric_limits<Real>::epsilon()) {
        std::cout << P << " and " << O << " already fall on a line in the direction of " << D << "\n";
        Real t = dot(PmO, D)/dot(D,D);
        if (t < bound[2].first || t > bound[2].second) {
            return vec<Real>(special_vec::NaNs);
        }
        return vec<Real>(u0, v0, t);
    }
    std::cout << "Plane normal before normalization = " << plane_normal << "\n";
    normalize(plane_normal);

    std::cout << "Our camera is at origin " << O << " and points in direction " << D << "\n";
    std::cout << "The plane normal is " << plane_normal << "\n";

    vec<Real> sigma_u = J.column(0);
    vec<Real> sigma_v = J.column(1);
    vec<Real> surface_normal = cross(sigma_u, sigma_v);
    normalize(surface_normal);
    std::cout << "The surface normal is " << surface_normal << "\n";
    

    vec<Real> tangent = cross(plane_normal, surface_normal);
    std::cout << "Tangent, prenormalization = " << tangent << "\n";
    normalize(tangent);
    std::cout << "Tangent, postnormalization = " << tangent << "\n";

    vec<Real,2> b;
    b[0] = dot(tangent, sigma_u);
    b[1] = dot(tangent, sigma_v);
    matrix<Real,2,2> M;
    M(0,0) = dot(sigma_u, sigma_u);
    M(0,1) = dot(sigma_u, sigma_v);
    M(1,0) = M(0,1);
    M(1,1) = dot(sigma_v, sigma_v);
    vec<Real,2> w = M.solve(b);
    Real uprime = w[0];
    Real vprime = w[1];

    std::cout << "u' = " << uprime << ", v' = " << vprime << "\n";

    // Second fundamental forms:
    Real e = 0;
    for (int64_t i = 0; i < 3; ++i) {
        e += H(i,0,0)*surface_normal[i];
    }

    Real f = 0;
    for (int64_t i = 0; i < 3; ++i) {
        f += H(i,0,1)*surface_normal[i];
    }

    Real g = 0;
    for (int64_t i = 0; i < 3; ++i) {
        g += H(i,1,1)*surface_normal[i];
    }

    Real kappa_b = e*uprime*uprime + 2*f*uprime*vprime + g*vprime*vprime;
    Real costheta = dot(surface_normal, plane_normal);
    if (costheta >= 1) {
        std::cerr << __FILE__ << ":" << __LINE__ << " The plane normal and the surface normal are parallel.\n";
        std::cerr << "You need to implement this degenerate case.\n";
    }
    Real s_ = kappa_b/(1-costheta*costheta);
    vec<Real> kappaN = s_*(-costheta*plane_normal + surface_normal);
    Real kappa = norm(kappaN);
    vec<Real> n = kappaN/kappa;


    Real a = kappa/2;
    Real tmp = dot(D,n);
    Real dt = dot(D,tangent);
    Real b_ = -tmp/dt;
    Real c = tmp*dot(O-P, tangent)/dt - dot(O-P,n);

    auto roots = quadratic_roots(a,b_,c);

    if (roots.size() == 0) {
        std::cerr << "No real roots found.\n";
    }
    else {
        std::cout << "Roots are " << roots[0] << ", " << roots[1] << "\n";
    }

    Real s = roots[0];
    bool unacceptable = false;
    if (u0 + s*uprime < bound[0].first || u0 + s*uprime > bound[0].second) {
        unacceptable = true;
        s = roots[1];
        if (u0 + s*uprime < bound[0].first || u0 + s*uprime > bound[0].second) {
            // Thrown out of bounds:
            return vec<Real>(special_vec::NaNs);
        }
    }

    if (v0 + s*vprime < bound[1].first || v0 + s*vprime > bound[1].second) {
        if (unacceptable) {
            return vec<Real>(special_vec::NaNs);
        }
        s = roots[1];
        if (v0 + s*vprime < bound[1].first || v0 + s*vprime > bound[1].second) {
            // Thrown out of bounds:
            return vec<Real>(special_vec::NaNs);
        }
    }
    u0 += s*uprime;
    v0 += s*vprime;

    std::cout << "New guess is (u0,v0) = (" << u0 << ", " << v0 << ")\n";


    P = h(u0,v0, J, H);
    std::cout << "New intersection point is " << P << "\n";
    // Now choose among the roots to find minimal t solution:

    return vec<Real>(special_vec::NaNs);
}

}
#endif