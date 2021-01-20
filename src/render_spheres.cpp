#include <drt/vec.hpp>

int main() {
    using Real = double;

    drt::vec<Real, 3> v(1.2, 3.4, 7.5);

    std::cout << v << "\n";
}