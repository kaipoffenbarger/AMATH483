#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

void dgemv(double a, const std::vector<std::vector<double>> &A, std::vector<double> &x, double b, std::vector<double> &y) {

    //checks
    if (A.empty() || A[0].empty()) {
        throw std::invalid_argument("Matrix A must not be empty.");
    }

    if (x.size() != A[0].size()) {
        throw std::invalid_argument("Size of x must equal the number of columns in A.");
    }

    if (y.size() != A.size()) {
        throw std::invalid_argument("Size of y must equal the number of rows in A.");
    }

    //calculation
    for (std::size_t i = 0; i < A.size(); ++i){
        double Ax_product_value = 0.0;
        
        for (std::size_t j = 0; j < A[i].size(); ++j){
            Ax_product_value += A[i][j] * x[j];
        }

        y[i] = a * Ax_product_value + b * y[i];
    }
}