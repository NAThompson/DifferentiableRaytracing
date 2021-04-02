#ifndef DRT_KO_METHOD_HPP
#define DRT_KO_METHOD_HPP
#include <drt/vec.hpp>
#include <drt/ray.hpp>
#include <drt/hittable.hpp>

namespace drt {

template<typename Real>
vec<Real,3> ko_method(hittable<Real>& h, ray<Real> const & r, Real u0, Real v0)
{
    vec<Real> O = r.origin();
    vec<Real> D = r.direction();
    matrix<Real,3,3> J;
    tensor<Real, 3, 2, 2> H;
    vec<Real> P = h(u0,v0, J, H);

    vec<Real> plane_normal = cross((P - O), D);
    normalize(plane_normal);

    vec<Real> sigma_u = J.column(0);
    vec<Real> sigma_v = J.column(1);
    vec<Real> surface_normal = cross(sigma_u, sigma_v);
    normalize(surface_normal);
    

    vec<Real> tangent = cross(plane_normal, surface_normal);
    normalize(tangent);

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
    Real s = kappa_b/(1-costheta*costheta);
    vec<Real> kappaN = s*(-costheta*plane_normal + surface_normal);
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

    // Now choose among the roots to find minimal t solution:

    return vec<Real>(special_vec::NaNs);
}

}
#endif