#ifndef DRT_HELICOID_HPP
#define DRT_HELICOID_HPP

#include <drt/hittable.hpp>
#include <drt/vec.hpp>
#include <drt/matrix.hpp>
#include <drt/roots.hpp>
#include <drt/newton.hpp>
#include <drt/halley.hpp>
#include <drt/tensor.hpp>
#include <drt/bounds.hpp>


namespace drt {

namespace {
    static int64_t enneper_error_count = 0;
    static int64_t enneper_hits = 0;
    static int64_t enneper_misses = 0;
}

// Parametrized via σ(u,v) = (u-u³/3 + uv², v-v³/3 + vu², u² - v²)
// Mathematica: ParametricPlot3D[{u - u^3/3 + u*v^2, v - v^3/3 + v*u^2, u^2 - v^2}, {u, -4, 4}, {v, -4, 4}]
template<typename Real>
class enneper : public hittable<Real> {
public:

    // We'll bound the enneper surface with a sphere of radius R just to make things sane.
    // The actual surface is not compact.
    enneper(Real radius)
       : radius_(radius)
    {
        if (radius_ <= 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << " Radius > 0 is required for a enneper.\n";
        }
    };

    virtual bool hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const override;

    virtual bool bounding_box(aabb<Real>& output_box) const override;

    virtual ~enneper() {
        #ifdef DEBUG
        std::cerr << "Enneper error count = " << enneper_error_count << "\n";
        std::cerr << "Enneper hits = " << enneper_hits << "\n";
        std::cerr << "Enneper misses = " << enneper_misses << "\n";
        #endif
    }

    inline vec<Real,3> operator()(Real u, Real v) const override {
        vec<Real, 3> p;
        p[0] = u*(1-u*u/3 + v*v);
        p[1] = v*(1-v*v/3 + u*u);
        p[2] = u*u - v*v;
        return p;
    }

    inline vec<Real,3> operator()(Real u, Real v, matrix<Real,3,3>& J) const {
        vec<Real> p = (*this)(u,v);
        // Note that a function σ:ℝ² -> ℝ³, the Jacobian is a 3x2 matrix.
        // However, that's pretty much useless for the purpose people actually use Jacobians for.
        // Hence, this returns a 3x3 matrix. The final column is filled with nans,
        // and must be populated by the user. The most common filling would be with -d for a ray
        // with r(t) = o + td.
        J(0,2) = std::numeric_limits<Real>::quiet_NaN();
        J(1,2) = std::numeric_limits<Real>::quiet_NaN();
        J(2,2) = std::numeric_limits<Real>::quiet_NaN();

        J(0,0) = 1 - p[2];
        J(1,0) = 2*u*v;
        J(2,0) = 2*u;

        J(0,1) = 2*u*v;
        J(1,1) = 1 + p[2];
        J(2,1) = -2*v;

        return p;
    }

    // I had an if constexpr that removed this code duplication.
    // I'm very displeased it played poorly with inheritance.
    inline vec<Real,3> operator()(Real u, Real v, matrix<Real,3,3>& J, tensor<Real, 3, 2, 2>& H) const {
        vec<Real> p = (*this)(u,v,J);
        // H(i,j,k) := ∂ⱼ∂ₖσᵢ = H(i,k,j).
        H(0,0,0) = -2*u;
        H(1,0,0) = 2*v;
        H(2,0,0) = 2;

        H(0,0,1) = 2*v;
        H(0,1,0) = H(0,0,1);

        H(1,0,1) = 2*u;
        H(1,1,0) = H(1,0,1);

        H(2,0,1) = 0;
        H(2,1,0) = 0;

        H(0,1,1) = 2*u;
        H(1,1,1) = -2*v;
        H(2,1,1) = -2;
        return p;
    }

    std::pair<Real,Real> gaussian_curvature_bounds() const {
        // K = -4/(1+u^2 + v^2)^4. u \in (-2,2), v \in (-2,2).
        return std::make_pair<Real, Real>(-4, -Real(4)/std::pow(Real(9),4));
    }


private:
    // Only used for the spherical bounding volume:
    Real radius_;

    vec<Real,3> compute_uvt(vec<Real> const & o, vec<Real> const & d, Real t_min, Real t_max) const
    {
        std::function<std::pair<vec<Real, 3>, matrix<Real,3,3>>(vec<Real, 3>)> f = [&o, &d, this](vec<Real, 3> w) {
            // w = (u,v,t).
            matrix<Real,3,3> J;
            vec<Real, 3> p = (*this)(w[0], w[1], J);

            J(0,2) = -d[0];
            J(1,2) = -d[1];
            J(2,2) = -d[2];
            return std::make_pair(p - o - w[2]*d, J);
        };

        std::array<vec<Real,3>, 12> sols;
        bounds<Real,3> bound({-2,0}, {-2,0}, {t_min, t_max});
        sols[0] = newton(f, bound, bound.center());
        bound[0] = std::make_pair<Real,Real>(0,2);
        bound[1] = std::make_pair<Real,Real>(0,2);
        sols[1] = newton(f, bound, bound.center());
        bound[0] = std::make_pair<Real,Real>(-1,1);
        bound[1] = std::make_pair<Real,Real>(-1,1);
        sols[2] = newton(f, bound, bound.center());
        bound[0] = std::make_pair<Real,Real>(-2,2);
        bound[1] = std::make_pair<Real,Real>(-2,2);

        for (size_t i = 3; i < sols.size(); ++i) {
            sols[i] = newton(f, bound, bound.random());
        }

        size_t idx = 0;
        Real tmin = std::numeric_limits<Real>::max();
        for (size_t i = 0; i < sols.size(); ++i) {
            vec<Real,3> sol = sols[i];
            if (sol[2] < tmin) {
                idx = i;
                tmin = sol[2];
            }
        }

        return sols[idx];
    }
};

template<typename Real>
bool enneper<Real>::hit(const ray<Real>& r, Real t_min, Real t_max, hit_record<Real>& rec) const {
    using std::abs;
    vec<Real, 3> oc = r.origin();
    Real a = drt::squared_norm(r.direction());
    Real b = 2*drt::dot(oc, r.direction());
    Real c = drt::squared_norm(oc) - radius_*radius_;
    auto roots = quadratic_roots(a, b, c);
    if (roots.size() != 2) {
        return false;
    }
    if (roots[1] < t_min) {
        return false;
    }
    if (roots[0] > t_max) {
        return false;
    }

    t_min = std::max(t_min, roots[0]);
    t_max = std::min(t_max, roots[1]);
    auto w = compute_uvt(r.origin(), r.direction(), t_min, t_max);
    if (std::isnan(w[0])) {
        ++enneper_misses;
        return false;
    }
    Real u = w[0];
    Real v = w[1];
    Real t = w[2];
    if(t < t_min || t > t_max) {
        ++enneper_misses;
        return false;
    }
    rec.t = t;
    rec.p = r(t);
    rec.u = u;
    rec.v = v;

    auto p = (*this)(u, v);

    Real residual = norm(rec.p - p);
    if (residual > 0.00001) {
        std::cerr << "Residual is " << residual << "\n";
        ++enneper_misses;
        ++enneper_error_count;
        return false;
    }

    vec<Real> dsigmadu(1-u*u + v*v, 2*u*v, 2*u);
    vec<Real> dsigmadv(2*u*v, 1-v*v + u*u, -2*v);
 
    vec<Real> gradient = drt::cross(dsigmadu, dsigmadv);
    rec.gradient_magnitude = norm(gradient);
    vec<Real> outward_normal = gradient/rec.gradient_magnitude;
    rec.set_face_normal(r, outward_normal);
    // E = (1+u^2+v2)^2.
    rec.E = squared_norm(dsigmadu);
    rec.F = 0;
    // G = (1+u^2+v2)^2.
    rec.G = squared_norm(dsigmadv);
    // Note the difference between this and https://mathworld.wolfram.com/EnnepersMinimalSurface.html
    // The sign of y is different in this convention.
    rec.e = 2;
    rec.f = 0;
    rec.g = -2;

    ++enneper_hits;
    return true;
}

template<typename Real>
bool enneper<Real>::bounding_box(aabb<Real>& output_box) const {
    output_box = aabb<Real>(
        vec<Real>(-radius_, -radius_, -radius_),
        vec<Real>(radius_, radius_, radius_));
    return true;
}

}
#endif
