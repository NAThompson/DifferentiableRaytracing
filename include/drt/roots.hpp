#ifndef DRT_ROOTS_HPP
#define DRT_ROOTS_HPP
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <optional>

namespace drt {

// https://stackoverflow.com/questions/48979861/numerically-stable-method-for-solving-quadratic-equations/50065711
template<class Real>
inline Real discriminant(Real const a, Real const b, Real const c)
{
   Real w = 4 * a * c;
   Real e = std::fma(-c, 4 * a, w);
   Real f = std::fma(b, b, -w);
   return f + e;
}

// Solves ax² + bx + c = 0.
// Only returns the real roots.
template<typename Real>
std::vector<Real> quadratic_roots(Real const a, Real const b, Real const c)
{
    using std::copysign;
    using std::sqrt;
    using std::abs;
    std::vector<Real> roots(2, std::numeric_limits<Real>::quiet_NaN());
    if (abs(a) <= std::numeric_limits<Real>::min())
    {
        if (b == 0 && c != 0)
        {
            if (c != 0)
            {
                roots.resize(0);
                return roots;
            }
            else
            {
                // Technically there are infinitely many roots.
                // Hopefully this interacts gracefully with intended use.
                roots[0] = 0;
                roots[1] = 0;
                return roots;
            }
        }
        // Degenerate roots:
        roots[0] = -c/b;
        roots[1] = -c/b;
        return roots;
    }
    if (abs(b) <= std::numeric_limits<Real>::min())
    {
        Real x0_sq = -c / a;
        if (x0_sq < 0)
        {
            roots.resize(0);
            return roots;
        }
        Real x0 = sqrt(x0_sq);
        roots[0] = -x0;
        roots[1] = x0;
        return roots;
    }
    Real delta = discriminant(a, b, c);
    if (delta < 0)
    {
        // All roots are complex:
        roots.resize(0);
        return roots;
    }
    Real q = -(b + copysign(sqrt(delta), b)) / 2;
    Real x0 = q / a;
    Real x1 = c / q;
    if (x0 < x1)
    {
        roots[0] = x0;
        roots[1] = x1;
        return roots;
    }
    
    roots[0] = x1;
    roots[1] = x0;
    return roots;
}

// Is a std::optional really better than a nan? I'm not so sure . . .
template<typename Real>
std::optional<Real> first_quadratic_root_in_range(Real a, Real b, Real c, Real tmin, Real tmax) {
    auto roots = quadratic_roots(a, b, c);
    if (roots.size() != 2) {
        return std::nullopt;
    }
    if (roots[0] >= tmin && roots[0] <= tmax) {
        return roots[0];
    }
    if (roots[1] >= tmin && roots[1] <= tmax) {
        return roots[1];
    }
    return std::nullopt;
}

template <typename Real> int sgn(Real val) {
    return (Real(0) < val) - (val < Real(0));
}

// Solves ax³ + bx² + cx + d = 0.
// Only returns the real roots, as these are the only roots of interest in ray intersection problems.
// Follows Numerical Recipes, Chapter 5, section 6.
template<typename Real>
std::vector<Real> cubic_roots(Real a, Real b, Real c, Real d) {
    using std::sqrt;
    using std::acos;
    using std::cos;
    using std::cbrt;
    using std::abs;
    if (abs(a) <= std::numeric_limits<Real>::min()) {
        return quadratic_roots(b, c, d);
    }
    if (abs(d) <= std::numeric_limits<Real>::min()) {
        auto v = quadratic_roots(a, b, c);
        v.push_back(0);
        std::sort(v.begin(), v.end());
        return v;
    }
    std::vector<Real> roots(3, std::numeric_limits<Real>::quiet_NaN());
    Real p = b/a;
    Real q = c/a;
    Real r = d/a;
    Real Q = (p*p - 3*q)/9;
    Real R = (2*p*p*p - 9*p*q + 27*r)/54;
    if (R*R < Q*Q*Q) {
        Real rtQ = sqrt(Q);
        Real theta = acos(R/(Q*rtQ));
        roots[0] = -2*rtQ*cos(theta/3) - p/3;
        roots[1] = -2*rtQ*cos((theta + 2*M_PI)/3) - p/3;
        roots[2] = -2*rtQ*cos((theta - 2*M_PI)/3) - p/3;
        std::sort(roots.begin(), roots.end());
        return roots;
    }
    // There is only 1 real root:
    roots.resize(1);
    Real arg = R*R - Q*Q*Q;
    Real A = -sgn(R)*cbrt(abs(R) + sqrt(arg));
    if (abs(A) <= std::numeric_limits<Real>::min()) {
        roots[0] = A - p/3;
    } else {
        roots[0] = A + Q/A - p/3;
    }
    return roots;
}


// Solves ax⁴ + bx³ + cx² + dx + e = 0.
// Only returns the real roots, as these are the only roots of interest in ray intersection problems.
// Follows Graphics Gems V: https://github.com/erich666/GraphicsGems/blob/master/gems/Roots3And4.c
template<typename Real>
std::vector<Real> quartic_roots(Real a, Real b, Real c, Real d, Real e) {
    using std::abs;
    using std::sqrt;
    if (std::abs(a) <= std::numeric_limits<Real>::min()) {
        return cubic_roots(b, c, d, e);
    }
    if (std::abs(e) <= std::numeric_limits<Real>::min()) {
        auto v = cubic_roots(a, b, c, d);
        v.push_back(Real(0));
        return v;
    }
    // Now solve x⁴ + Ax³ + Bx² + Cx + D = 0.
    Real A = b/a;
    Real B = c/a;
    Real C = d/a;
    Real D = e/a;
    Real Asq = A*A;
    // Let x = y - A/4:
    // Mathematica: Expand[(y - A/4)^4 + A*(y - A/4)^3 + B*(y - A/4)^2 + C*(y - A/4) + D]
    // We now solve y⁴ + py² + qy + r = 0.
    Real p = B - 3*Asq/8;
    Real q = C - A*B/2 + Asq*A/8;
    Real r = D - A*C/4 + Asq*B/16 - 3*Asq*Asq/256;
    if (std::abs(r) <= std::numeric_limits<Real>::min()) {
        auto v = cubic_roots(Real(1), Real(0), p, q);
        v.push_back(Real(0));
        for (auto & y : v) {
            y -= A/4;
        }
        std::sort(v.begin(), v.end());
        return v;
    }
    // Biquadratic case:
    if (std::abs(q) <= std::numeric_limits<Real>::min()) {
        auto w = quadratic_roots(Real(1), p, r);
        std::vector<Real> v;
        for (auto r : w) {
            if (r >= 0) {
                Real rtr = sqrt(r);
                v.push_back(rtr - A/4);
                v.push_back(-rtr - A/4);
            }
        }
        std::sort(v.begin(), v.end());
        return v;
    }

    // Now split the depressed cubic into two quadratics:
    // y⁴ + py² + qy + r = (y² + sy + u)(y² - sy + v) = y⁴ + (v+u-s²)y² + s(v - u)y + uv
    // So p = v+u-s², q = s(v - u), r = uv.
    // Then (v+u)² - (v-u)² = 4uv = 4r = (p+s²)² - q²/s².
    // Multiply through by s² to get s²(p+s²)² - q² - 4rs² = 0, which is a cubic in s².
    // Then we let z = s², to get
    // z³ + 2pz² + (p² - 4r)z - q² = 0.
    auto z_roots = cubic_roots(Real(1), 2*p, p*p - 4*r, -q*q);
    // z = s², so s = sqrt(z).
    // No real roots:
    if (z_roots.back() <= 0) {
        std::vector<Real> v(0);
        return v;
    }
    Real s = std::sqrt(z_roots.back());

    // s is nonzero, because we took care of the biquadratic case.
    Real v = (p + s*s + q/s)/2;
    Real u = v - q/s;
    // Now solve y² + sy + u = 0:
    auto roots1 = quadratic_roots(Real(1), s, u);

    // Now solve y² - sy + v = 0:
    auto roots2 = quadratic_roots(Real(1), -s, v);
    for (auto root : roots2) {
        roots1.push_back(root);
    }

    for (auto& r : roots1) {
        r -= A/4;
    }

    // This is not super accurate. Clean up the roots with a Halley iterate.
    for (auto &r : roots1) {
        Real df = 4*a*r + 3*b;
        df = df*r + 2*c;
        df = df*r + d;
        Real d2f = 12*a*r + 6*b;
        d2f = d2f*r + 2*c;
        Real f = a*r + b;
        f = f*r + c;
        f = f*r + d;
        f = f*r + e;
        Real denom = 2*df*df - f*d2f;
        if (std::abs(denom) > std::numeric_limits<Real>::min())
        {
            r -= 2*f*df/denom;
        }
    }

    std::sort(roots1.begin(), roots1.end());
    return roots1;
}


template<typename Real>
std::pair<Real, Real> bisect(std::function<Real(Real)> f, Real tmin, Real tmax)
{
    if (tmax < tmin) {
        std::cerr << "Arguments to bisect in wrong order: [tmin, tmax] = [" << tmin << ", " << tmax << "]\n";
    }
    Real f_min = f(tmin);
    if (f_min == 0) {
        return std::make_pair(tmin, tmin);
    }

    // This is a pretty expensive check, but nonexistence of roots is important!
    Real dt = std::min((tmax - tmin)/8, Real(1));
    if (dt < 10*std::numeric_limits<Real>::epsilon()*abs(tmin)) {
        dt = 10*std::numeric_limits<Real>::epsilon()*abs(tmin);
    }
    Real f_max = f(tmin + dt);
    while (f_min*f_max > 0 && tmin < tmax) {
        f_min = f_max;
        tmin += dt;
        f_max = f(tmin + dt);
    }
    // No sign change found anywhere!
    if (tmin >= tmax)
    {
        return std::make_pair(std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN());
    }
    tmax = tmin + dt;

    int count = 0;
    while (count++ < 40) {
        Real mid = (tmin + tmax) / 2;
        Real fmid = f(mid);
        if ((mid == tmax) || (mid == tmin))
        {
            break;
        }

        if (fmid == 0)
        {
            tmin = tmax = mid;
            break;
        }
        if (f_min*f_min < 0)
        {
            tmax = mid;
        }
        else
        {
            tmin = mid;
            f_min = fmid;
        }
    }
    return std::make_pair(tmin, tmax);
}

}
#endif
