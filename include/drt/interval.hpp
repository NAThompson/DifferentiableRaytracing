#ifndef DRT_INTERVAL_HPP
#define DRT_INTERVAL_HPP
#include <iostream>
#include <cfenv>

namespace drt {

template<typename Real>
class interval {
public:
    interval() : a_{std::numeric_limits<Real>::quiet_NaN()}, b_{std::numeric_limits<Real>::quiet_NaN()}
    {}

    interval(Real a, Real b) : a_{a}, b_{b} {
        if (b_ < a_) {
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << " An empty interval has been created.\n";
            std::cerr << "a <= b is require for an interval [a,b], but obtained a = " << a << " and b = " << b << "\n";
        }
    }

    friend std::ostream& operator<<(std::ostream& out, interval<Real> const & i)
    {
        out << "[" << i.a_ << ", " << i.b_ << "]";
        return out;
    }

    inline bool empty() const {
        return b_ < a_;
    }

    inline Real lower() const {
        return a_;
    }

    inline Real upper() const {
        return b_;
    }

    inline Real width() const {
        return b_ - a_;
    }

    inline Real degenerate() const {
        return a_ == b_;
    }

    inline interval<Real>& operator-=(const interval<Real>& r)
    {
        const int originalRounding = fegetround();
        std::fesetround(FE_DOWNWARD);
        a_ -= r.b_;
        std::fesetround(FE_UPWARD);
        b_ -= r.a_;
        std::fesetround(originalRounding);
        return *this;
    }

    inline bool operator>=(const interval<Real>& r) const;

    inline interval<Real> operator+(const interval<Real>& r) const;

    Real a_;
    Real b_;
};


template<typename Real>
inline bool interval<Real>::operator>= (const interval<Real>& r) const
{
  if (!empty()) 
  {
    if (a_ >= r.b_)
    {
        return true;
    }
    else if (b_ < r.a_)
    {
        return false;
    }
  }
  std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << " Non-comparable intervals.\n";
  std::cerr << "\tExpected [" << a_ << ", " << b_ << "] >= ["  << r.a_ << ", " << r.b_ << "]\n";
  return false;
}

template<typename Real>
inline interval<Real> interval<Real>::operator+(const interval<Real>& r) const
{
    // Need to check the speed of resetting the rounding mode.
    const int originalRounding = fegetround();
    std::fesetround(FE_DOWNWARD);
    Real a = a_ + r.a_;
    std::fesetround(FE_UPWARD);
    Real b = b_ + r.b_;
    std::fesetround(originalRounding);
    return interval<Real>(a, b);
}

template<typename Real>
inline interval<Real> operator-(const interval<Real>& x, const interval<Real>& y)
{
    const int originalRounding = fegetround();
    std::fesetround(FE_DOWNWARD);
    Real a = x.lower() - y.upper();
    std::fesetround(FE_UPWARD);
    Real b = x.upper() - y.lower();
    std::fesetround(originalRounding);
    return interval<Real>(a, b);
}

template<typename Real>
inline interval<Real> operator*(const interval<Real>& x, const interval<Real>& y)
{
    std::array<Real, 4> arr{x.lower()*y.lower(), x.lower()*y.upper(), x.upper()*y.lower(), x.upper()*y.upper()};
    Real a = std::numeric_limits<Real>::infinity();
    Real b = -std::numeric_limits<Real>::infinity();

    for (size_t i = 0; i < 4; ++i) {
        if (a > arr[i]) {
            a = arr[i];
        }
        if (b < arr[i]) {
            b = arr[i];
        }
    }
    return interval<Real>(a, b);
}


template<typename Real>
inline interval<Real> operator/(const interval<Real>& x, const interval<Real>& y)
{
    if (y.lower() <= 0 && y.upper() > 0) {
        // Division by zero in intervals is just fine,
        // the problem is that they undergo mitosis.
        if (y.lower() < 0) {
            std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << " Division by an interval containing zero.\n";
            std::cerr << "\tYou need a mechanism to deal with multi-intervals.\n";
        }
        return x*interval<Real>(1/y.upper(), std::numeric_limits<Real>::infinity());
    }
    interval<Real> inv_y(1/y.upper(), 1/y.lower());
    return x*inv_y;
}

template<typename Real>
inline interval<Real> pow(const interval<Real>& x, int p)
{
    if ( ((p & 1) == 0) && x.lower() <= 0)
    {
        return interval<Real>(0, std::pow(x.upper(), p));
    }
    return interval<Real>(std::pow(x.lower(), p), std::pow(x.upper(), p));
}

template<typename Real>
inline interval<Real> abs(const interval<Real>& x)
{
    Real a = min(abs(x.lower()), abs(x.upper()));
    if (x.lower() < 0 && x.upper() > 0) {
        a = 0;
    }
    Real b = max(abs(x.lower()), abs(x.upper()));
    return interval<Real>(a,b);
}

template<typename Real>
inline interval<Real> intersection(const interval<Real>& x, const interval<Real>& y)
{
    if ( (y.lower() > x.upper()) || (x.lower() > y.upper())) {
        // Default contstructor is nan's; is that a reasonable empty set?
        return interval<Real>();
    }
    Real a = std::max(x.lower(), y.lower());
    Real b = std::min(x.upper(), y.upper());
    return interval<Real>(a,b);
}

// The union of two intervals is no necessarily an interval; there could be a gap in the middle.
// The hull is the union + gaps.
template<typename Real>
inline interval<Real> hull(const interval<Real>& x, const interval<Real>& y)
{
    Real a = std::min(x.lower(), y.lower());
    Real b = std::min(x.upper(), y.upper());
    return interval<Real>(a,b);
}

template<typename Real>
inline interval<Real> exp(const interval<Real>& x)
{
    const int originalRounding = fegetround();
    std::fesetround(FE_DOWNWARD);
    Real a = std::exp(x.lower());
    std::fesetround(FE_UPWARD);
    Real b = std::exp(x.upper());
    std::fesetround(originalRounding);
    return interval<Real>(a,b);
}

template<typename Real>
inline interval<Real> log(const interval<Real>& x)
{
    const int originalRounding = fegetround();
    std::fesetround(FE_DOWNWARD);
    Real a = std::log(x.lower());
    std::fesetround(FE_UPWARD);
    Real b = std::log(x.upper());
    std::fesetround(originalRounding);
    return interval<Real>(a,b);
}

template<typename Real>
inline interval<Real> sqrt(const interval<Real>& x)
{
    const int originalRounding = fegetround();
    std::fesetround(FE_DOWNWARD);
    // If x.lower() < 0, this spits out a NaN.
    // There are policies for this that some libraries set out,
    // for example, you can set the a = 0 if x.lower() < 0.
    // This makes some computations proceed without harm,
    // but also makes others wrong!
    Real a = std::sqrt(x.lower());
    std::fesetround(FE_UPWARD);
    Real b = std::sqrt(x.upper());
    std::fesetround(originalRounding);
    return interval<Real>(a,b);
}

}
#endif
