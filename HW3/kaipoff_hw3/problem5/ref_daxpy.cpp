#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

void daxpy(double a, const std::vector<double> &x, std::vector<double> &y) {

    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors are different lengths! Vectors must be the same size");
    }

    for (std::size_t i = 0; i < x.size(); ++i){
        y[i] = a * x[i] + y[i];
    }
}