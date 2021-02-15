#ifndef DRT_HITTABLE_HPP
#define DRT_HITTABLE_HPP
#include <limits>
#include <drt/roots.hpp>
#include <drt/ray.hpp>
#include <drt/vec.hpp>
#include <drt/aabb.hpp>


namespace drt {

template<typename Real>
class material;

template<typename Real>
struct hit_record {
    vec<Real, 3> p;
    vec<Real, 3> normal;
    Real gradient_magnitude = std::numeric_limits<Real>::quiet_NaN();
    Real t = std::numeric_limits<Real>::quiet_NaN();
    // Parametric coordinates on the surface.
    // Not used by everything, but if it is used, u, v \in [0,1].
    Real u = std::numeric_limits<Real>::quiet_NaN();
    Real v = std::numeric_limits<Real>::quiet_NaN();
    bool front_face;
    std::shared_ptr<material<Real>> mat_ptr;

    inline void set_face_normal(const ray<Real>& r, const vec<Real, 3>& outward_normal) {
        front_face = drt::dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }

    // First fundamental form. Computed under the assumption of [0,1] parametrization.
    Real E = std::numeric_limits<Real>::quiet_NaN();
    Real F = std::numeric_limits<Real>::quiet_NaN();
    Real G = std::numeric_limits<Real>::quiet_NaN();

    // Second fundamental form. Again, [0,1] parametrization assumed.
    Real L = std::numeric_limits<Real>::quiet_NaN();
    Real M = std::numeric_limits<Real>::quiet_NaN();
    Real N = std::numeric_limits<Real>::quiet_NaN();

    // Pressley, Elementary Differential Geometry, Corollary 8.1.3
    Real gaussian_curvature() const {
        return (L*N - M*M)/(E*G - F*F);
    }

    Real mean_curvature() const {
        return (L*G - 2*M*F + N*E)/(2*(E*G - F*F));
    }

    std::pair<Real, Real> principal_curvatures() const {
        Real H = mean_curvature();
        Real K = gaussian_curvature();
        auto v = quadratic_roots(Real(1), -2*H, K);
        return std::make_pair(v[0], v[1]);
    }

    friend std::ostream& operator<<(std::ostream & os, const hit_record<Real> & rec) {
        os << "Ray intersects ";
        if (rec.front_face) {
            os << "front face";
        } else {
            os << "back face";
        }
        os << " of object at " << rec.p << " and time " << rec.t << "\n";
        os << "Normal is " << rec.normal << " with magnitude of gradient = " << rec.gradient_magnitude << "\n";
        os << "Parametric coordinates: (u,v) = (" << rec.u << ", " << rec.v << ")\n";
        os << "First fundamental form (E,F,G) = (" << rec.E << ", " << rec.F << ", " << rec.G << ")\n";
        os << "Second fundamental form (L, M, N) = (" << rec.L << ", " << rec.M << ", " << rec.N << ")\n";
        os << "Gaussian curvature K = " << rec.gaussian_curvature() << ", mean curvature H = " << rec.mean_curvature() << ".\n";
        auto [k1, k2] = rec.principal_curvatures();
        os << "Principle curvatures are κ₁ = " << k1 << " and κ₂ = " << k2 << "\n";
        return os;
    }
};

template<typename Real>
class hittable {
public:
    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const = 0;

    virtual bool bounding_box(aabb<Real>& output_box) const = 0;
};

}
#endif