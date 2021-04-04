#ifndef DRT_KO_METHOD_HPP
#define DRT_KO_METHOD_HPP
#include <drt/vec.hpp>
#include <drt/ray.hpp>
#include <drt/hittable.hpp>

namespace drt {

#define DRT_DEBUG_KO_METHOD 1

template<typename Real>
vec<Real,3> ko_method(hittable<Real>& h, ray<Real> const & r, bounds<Real,3> const & bound, Real u0 = std::numeric_limits<Real>::quiet_NaN(), Real v0 = std::numeric_limits<Real>::quiet_NaN())
{
    vec<Real> O = r.origin();
    vec<Real> D = r.direction();
    #if DRT_DEBUG_KO_METHOD
    std::cout << "The ray is " << r << "\n";
    #endif
    if (std::isnan(u0)) {
        u0 = (bound[0].first + bound[0].second)/2;
    }
    if (std::isnan(v0)) {
        v0 = (bound[1].first + bound[1].second)/2;
    }

    matrix<Real,3,3> J;
    tensor<Real, 3, 2, 2> H;
    bool converged = false;
    do {
    vec<Real> P = h(u0,v0, J, H);

    #if DRT_DEBUG_KO_METHOD
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    std::cout << "\n\tGuess (u₀,v₀) = (" << u0 << ", " << v0 << ") gives "
              << "σ(u₀,v₀) = " << P << "\n";
    #endif

    vec<Real> PmO = P - O;
    vec<Real> plane_normal = cross(PmO, D);

    // TODO: What is the real condition that these two fall on a line using the rounding model of floating point arithmetic?
    if (squared_norm(plane_normal) < std::numeric_limits<Real>::epsilon()*std::numeric_limits<Real>::epsilon()) {
        #if DRT_DEBUG_KO_METHOD
        std::cout << "\t" << P << " and " << O << " fall on a line in the direction of " << D << "; the iteration has converged.\n";
        #endif
        Real t = dot(PmO, D)/dot(D,D);
        if (t < bound[2].first || t > bound[2].second) {
            return vec<Real>(special_vec::NaNs);
        }
        return vec<Real>(u0, v0, t);
    }
    //std::cout << "Plane normal before normalization = " << plane_normal << "\n";
    normalize(plane_normal);
    //std::cout << "The plane normal is " << plane_normal << "\n";

    vec<Real> sigma_u = J.column(0);
    vec<Real> sigma_v = J.column(1);
    vec<Real> surface_normal = cross(sigma_u, sigma_v);
    if (squared_norm(surface_normal) < std::numeric_limits<Real>::epsilon()*std::numeric_limits<Real>::epsilon()) {
        std::cerr << __FILE__ << ":" << __LINE__ << "\n\tAt (u₀,v₀) = (" << u0 << ", " << v0 << "), the surface patch σ is irregular; "
                  << "∂ᵤσ⨯∂ᵥσ = " << surface_normal << " cannot be sensibly normalized and hence the surface has an ill-defined normal.\n";
        return vec<Real>(special_vec::NaNs);
    }
    //std::cout << "\tThe surface normal is " << surface_normal << "\n";
    normalize(surface_normal);
    //std::cout << "\tThe surface normal is " << surface_normal << "\n";
    

    vec<Real> tangent = cross(plane_normal, surface_normal);
   // std::cout << "\tTangent, prenormalization = " << tangent << "\n";
    if (squared_norm(tangent) < std::numeric_limits<Real>::epsilon()*std::numeric_limits<Real>::epsilon()) {
        std::cerr << __FILE__ << ":" << __LINE__ << " the tangent is " << tangent << ", which is too small!\n";
        std::cerr << "You must reckon with this degenerate case.\n";
        return vec<Real>(special_vec::NaNs);
    }
    normalize(tangent);
    //std::cout << "\tTangent, postnormalization = " << tangent << "\n";

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
    Real f = 0;
    Real g = 0;
    for (int64_t i = 0; i < 3; ++i) {
        e += H(i,0,0)*surface_normal[i];
        f += H(i,0,1)*surface_normal[i];
        g += H(i,1,1)*surface_normal[i];
    }

    Real kappa_b = e*uprime*uprime + 2*f*uprime*vprime + g*vprime*vprime;
    Real costheta = dot(surface_normal, plane_normal);
    if (costheta >= 1) {
        std::cerr << __FILE__ << ":" << __LINE__ << " The plane normal and the surface normal are parallel.\n";
        std::cerr << "You need to implement this degenerate case.\n";
        return vec<Real>(special_vec::NaNs);
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
        #if DRT_DEBUG_KO_METHOD
        std::cout << "\tNo real roots found; giving up.\n";
        #endif
        return vec<Real>(special_vec::NaNs);
    }

    vec<Real> sol0(u0 + roots[0]*uprime, v0 + roots[0]*vprime, (roots[0] + dot(PmO, tangent))/dt);
    vec<Real> sol1(u0 + roots[1]*uprime, v0 + roots[1]*vprime, (roots[1] + dot(PmO, tangent))/dt);
    Real s = roots[0];
    bool contains0 = bound.contains(sol0);
    bool contains1 = bound.contains(sol1);
    if (!contains0 && !contains1) {
        #if DRT_DEBUG_KO_METHOD
        std::cout << "\tNeither root generates an in-bounds update.\n";
        #endif
        return vec<Real>(special_vec::NaNs);
    }
    else if (!contains0) 
    {
        s = roots[1];
    }
    else if (!contains1) 
    {
        s = roots[0];
    }
    else {
        //
        if (sol1[2] < sol0[2]) {
            s = roots[1];
        }
    }

    converged = ((s*uprime)*(s*uprime) + (s*vprime)*(s*vprime)) < std::numeric_limits<Real>::epsilon()*std::numeric_limits<Real>::epsilon();
    u0 += s*uprime;
    v0 += s*vprime;

    #if DRT_DEBUG_KO_METHOD
    std::cout << "\tThe roots are s = " << roots[0] << " and " << roots[1] << ".\n";
    if (contains0 && contains1) {
        std::cout << "\tBoth roots generate in-bounds updates.\n";
    }
    else if (contains0)
    {
        std::cout << "\tOnly root " << roots[0] << " generates an in-bounds update.\n";
    }
    else if (contains1) {
        std::cout << "\tOnly root " << roots[1] << " generates an in-bounds update.\n";
    }
    std::cout << "\t(u̇,v̇) = (" << uprime << ", " << vprime << ") and s = " << s << ".\n";
    std::cout << "\tThis gives the correction (δu, δv) = (" << s*uprime << ", " << s*vprime << ").\n";
    #endif


    } while (!converged);

    return vec<Real>(special_vec::NaNs);
}

}
#endif