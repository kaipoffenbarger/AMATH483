#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

void dgemm(double a, const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B, double b, std::vector<std::vector<double>> &C) {

    //dimensions
    //rows
    std::size_t Arows = A.size();
    std::size_t Brows = B.size();
    std::size_t Crows = C.size();
    //cols
    std::size_t Acols = A[0].size();
    std::size_t Bcols = B[0].size();
    std::size_t Ccols = C[0].size();

    //checks
    if (A.empty() || A[0].empty()) {
        throw std::invalid_argument("Matrix A must not be empty.");
    }
    if (B.empty() || B[0].empty()) {
        throw std::invalid_argument("Matrix A must not be empty.");
    }
    if (C.empty() || C[0].empty()) {
        throw std::invalid_argument("Matrix A must not be empty.");
    }

    if (Acols != Brows) {
        throw std::invalid_argument("Column dimension of A must equal row dimension of B.");
    }
    if (Arows != Crows) {
        throw std::invalid_argument("Row dimension of A must equal row dimension of C.");
    }
    if (Bcols != Ccols) {
        throw std::invalid_argument("Column dimension of B must equal column dimension of C.");
    }

    //calculation
    for (std::size_t i = 0; i < Arows; ++i) {
        for (std::size_t j = 0; j < Bcols; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < Acols; ++k) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = a * sum + b * C[i][j];
        }
    }

}