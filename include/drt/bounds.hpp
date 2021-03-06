#ifndef DRT_BOUNDS_HPP
#define DRT_BOUNDS_HPP
#include <utility>
#include <drt/vec.hpp>

namespace drt {

template<typename Real, int64_t dimension>
class bounds {
public:
    bounds(std::pair<Real, Real> const & b)
    {
        static_assert(dimension == 1);
        bounds_[0] = b;
        assert(b.first <= b.second);
    }

    bounds(std::pair<Real, Real> const & b0, std::pair<Real, Real> const & b1)
    {
        static_assert(dimension == 2);
        bounds_[0] = b0;
        assert(b0.first <= b0.second);
        bounds_[1] = b1;
        assert(b1.first <= b1.second);
    }

    bounds(std::pair<Real, Real> const & b0, std::pair<Real, Real> const & b1, std::pair<Real, Real> const & b2)
    {
        static_assert(dimension == 3);
        bounds_[0] = b0;
        assert(b0.first <= b0.second);
        bounds_[1] = b1;
        assert(b1.first <= b1.second);
        bounds_[2] = b2;
        assert(b2.first <= b2.second);
    }

    inline bool contains(vec<Real, dimension> const & v) const
    {
        for (int64_t i = 0; i < dimension; ++i) {
            if (v[i] < bounds_[i].first || v[i] > bounds_[i].second) {
                return false;
            }
        }
        return true;
    }

    inline bool contains(Real x, Real y) const
    {
        static_assert(dimension == 2);
        if (x < bounds_[0].first || x > bounds_[0].second) {
            return false;
        }
        if (y < bounds_[1].first || x > bounds_[1].second) {
            return false;
        }
        return true;
    }

    inline bool contains(Real x, Real y, Real z) const
    {
        static_assert(dimension == 3);
        if (x < bounds_[0].first || x > bounds_[0].second) {
            return false;
        }
        if (y < bounds_[1].first || x > bounds_[1].second) {
            return false;
        }
        if (z < bounds_[2].first || z > bounds_[2].second) {
            return false;
        }
        return true;
    }

    vec<Real, dimension> center() const {
        vec<Real,dimension> v;
        for (int64_t i = 0; i < dimension; ++i) {
            v[i] = (bounds_[i].first + bounds_[i].second)/2;
        }
        return v;
    }

    inline Real volume() const {
        Real v = 1;
        for (int64_t i = 0; i < dimension; ++i) {
            v *= (bounds_[i].second - bounds_[i].first);
        }
        return v;
    }

    inline bool empty() const {
        return (this->volume() <= 0);
    }

    // This isn't fast and sucks you entropy pool dry!!!!
    vec<Real, dimension> random() const {
        std::random_device rd;
        vec<Real, dimension> v;
        for (int64_t i = 0; i < dimension; ++i) {
            std::uniform_real_distribution<Real> dis(bounds_[i].first, bounds_[i].second);
            v[i] = dis(rd);
        }
        return v;
    }

    friend std::ostream& operator<<(std::ostream & os, const bounds<Real, dimension> & b)
    {
        os << "{";
        for (int64_t i = 0; i < dimension - 1; ++i)
        {
            os << "[" << b.bounds_[i].first << ", " << b.bounds_[i].second << "], ";
        }
        os << "[" << b.bounds_.back().first << ", " << b.bounds_.back().second << "]}";
        return os;
    }

private:
    std::array<std::pair<Real, Real>, dimension> bounds_;
};

}
#endif