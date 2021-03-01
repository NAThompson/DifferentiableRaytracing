#ifndef DRT_TENSOR_HPP
#define DRT_TENSOR_HPP
#include <cmath>
#include <array>
#include <iostream>
#include <initializer_list>

namespace drt {

template<typename Real, int64_t rows = 3, int64_t cols = 3, int64_t panels = 3>
class tensor {
public:

    tensor() {};

    Real operator()(int64_t i, int64_t j, int64_t k) const {
        return T_[i*panels + j*cols + k];
    }

    Real& operator()(int64_t i, int64_t j, int64_t k) {
        return T_[i*panels + j*cols + k];
    }

    friend std::ostream& operator<<(std::ostream & os, const tensor<Real, rows, cols, panels> & v)
    {
        // This is totally wrong.
        for (int64_t i = 0; i < panels; ++i) {
            os << "[";
            for (int64_t j = 0; j < cols - 1; ++j) {
                os << v.T_[i*panels + j*cols] << ", ";
            }
            os << v.T_[i*panels + cols] << "]\n";
        }
        return os;
    }

private:
    Real T_[rows*cols*panels];
};


}
#endif
