#ifndef TENSOR_TEST_HPP
#define TENSOR_TEST_HPP
#include <drt/tensor.hpp>

using namespace drt;

TEST(TensorTest, TwoByTwo) {
    using Real = double;
    tensor<Real, 2,2,2> t;
    for (int64_t i = 0; i < 2; ++i) {
        for (int64_t j = 0; j < 2; ++j) {
            for (int64_t k = 0; k < 2; ++k) {
                t(i,j,k) = i*i + j*j + k*k;
            }
        }
    }
    std::cout << t << "\n";
}


TEST(TensorTest, ThreeByThree) {
    using Real = double;
    tensor<Real, 3,3,3> t;
    
}
#endif