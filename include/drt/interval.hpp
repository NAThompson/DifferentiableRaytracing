#ifndef DRT_INTERVAL_HPP
#define DRT_INTERVAL_HPP
#include <iostream>

namespace drt {

template<typename Real>
class interval {
public:
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
        return !(a_ <= b_);
    }

    inline Real lower() const {
        return a_;
    }

    inline Real upper() const {
        return b_;
    }

    inline bool operator>=(const interval<Real>& r) const;

    inline interval<Real> operator+(const interval<Real>& r) const;

private:
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
    return interval<Real>(a_ + r.a_, b_ + r.b_);
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
    if (y.lower() < 0 && y.upper() > 0) {
        std::cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ << " Division by an interval containing zero.\n";
        std::cerr << "\tYou need a mechanism to deal with multi-intervals.\n";
        return x*interval<Real>(1/y.upper(), std::numeric_limits<Real>::infinity());
    }
    interval<Real> inv_y(1/y.upper(), 1/y.lower());
    return x*inv_y;
}


}
#endif