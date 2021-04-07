#ifndef DRT_HITTABLE_HPP
#define DRT_HITTABLE_HPP
#include <limits>
#include <drt/roots.hpp>
#include <drt/ray.hpp>
#include <drt/vec.hpp>
#include <drt/matrix.hpp>
#include <drt/tensor.hpp>
#include <drt/aabb.hpp>


namespace drt {

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

    inline void set_face_normal(const ray<Real>& r, const vec<Real, 3>& outward_normal) {
        front_face = drt::dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }

    // First fundamental form. Computed under the assumption of [0,1] parametrization.
    // Let p = σ(u,v). Then
    // E = ‖∂σ/∂u‖², F = ∂σ/∂u·∂σ/∂v, G = ‖∂σ/∂v‖²
    Real E = std::numeric_limits<Real>::quiet_NaN();
    Real F = std::numeric_limits<Real>::quiet_NaN();
    Real G = std::numeric_limits<Real>::quiet_NaN();

    // Second fundamental form. Again, [0,1] parametrization assumed.
    // Follow the notation here: http://www.pbr-book.org/3ed-2018/Shapes/Spheres.html#PartialDerivativesofNormalVectors
    // e = ∂²σ/∂u²·n, f = ∂²σ/∂u∂v·n, g = ∂²σ/∂v²·n
    Real e = std::numeric_limits<Real>::quiet_NaN();
    Real f = std::numeric_limits<Real>::quiet_NaN();
    Real g = std::numeric_limits<Real>::quiet_NaN();

    Real condition_number = std::numeric_limits<Real>::quiet_NaN();
    // Pressley, Elementary Differential Geometry, Corollary 8.1.3
    Real gaussian_curvature() const {
        return (e*g - f*f)/(E*G - F*F);
    }

    Real mean_curvature() const {
        return (e*G - 2*f*F + g*E)/(2*(E*G - F*F));
    }

    std::pair<Real, Real> principal_curvatures() const {
        Real H = mean_curvature();
        Real K = gaussian_curvature();
        auto v = quadratic_roots(Real(1), -2*H, K);
        if (v.size() != 2) {
            return std::make_pair(std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN());
        }
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
        os << "Second fundamental form (e,f,g) = (" << rec.e << ", " << rec.f << ", " << rec.g << ")\n";
        os << "Gaussian curvature K = " << rec.gaussian_curvature() << ", mean curvature H = " << rec.mean_curvature() << ".\n";
        auto [k1, k2] = rec.principal_curvatures();
        os << "Principle curvatures are κ₁ = " << k1 << " and κ₂ = " << k2 << ".\n";
        return os;
    }
};

template<typename Real>
class hittable {
public:
    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const = 0;

    virtual bool bounding_box(aabb<Real>& output_box) const = 0;

    virtual vec<Real,3> operator()(Real u, Real v) const = 0;

    virtual vec<Real,3> operator()(Real, Real, matrix<Real,3,3> & J) const {
        for (int64_t i = 0; i < 3; ++i) {
            for (int64_t j = 0; j < 3; ++j) {
                J(i,j) = std::numeric_limits<Real>::quiet_NaN();
            }
        }
        return vec<Real>(special_vec::NaNs);
    };

    virtual vec<Real,3> operator()(Real, Real, matrix<Real,3,3> & J, tensor<Real,3,2,2> & H) const {
        for (int64_t i = 0; i < 3; ++i) {
            for (int64_t j = 0; j < 3; ++j) {
                J(i,j) = std::numeric_limits<Real>::quiet_NaN();
            }
        }
        for (int64_t i = 0; i < 3; ++i) {
            for (int64_t j = 0; j < 2; ++j) {
                for (int64_t k = 0; k < 2; ++k) {
                    H(i,j,k) = std::numeric_limits<Real>::quiet_NaN();
                }
            }
        }

        return vec<Real>(special_vec::NaNs);
    };

};

}
#endif