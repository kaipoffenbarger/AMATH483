#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

template <typename T>
void axpy(T a, const std::vector<T> &x, std::vector<T> &y) {

    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors are different lengths! Vectors must be the same size");
    }

    for (std::size_t i = 0; i < x.size(); ++i){
        y[i] = a * x[i] + y[i];
    }
}

template void axpy(double a, const std::vector<double> &x, std::vector<double> &y);
template void axpy(int a, const std::vector<int> &x, std::vector<int> &y);
template void axpy(float a, const std::vector<float> &x, std::vector<float> &y);