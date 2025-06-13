#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include "refBLAS.hpp"

template <typename T>
void displayMatrix(const std::vector<std::vector<T>>& matrix) {
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        for (std::size_t j =0; j < matrix[0].size(); ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {

    std::mt19937 rng(5); // Random seed
    std::uniform_int_distribution<int> dist(1, 10); // Non-zero values
    int n = 3;

    int alpha = dist(rng);
    int beta = dist(rng);

    std::vector<int> x(n), y(n);

    for (int i = 0; i < n; ++i) {
        x[i] = dist(rng);
        y[i] = dist(rng);
    }
            
    std::vector<std::vector<int>> A(n, std::vector<int>(n));
    std::vector<std::vector<int>> B(n, std::vector<int>(n));
    std::vector<std::vector<int>> C(n, std::vector<int>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = dist(rng);
            B[i][j] = dist(rng);
            C[i][j] = dist(rng);
        }
    }

    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "beta: " << beta << std::endl;
    
    std::cout << "Original y: ";

    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << y[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "x: ";

    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << y[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "A: ";
    displayMatrix(A);
    std::cout << std::endl;

    std::cout << "B: ";
    displayMatrix(B);
    std::cout << std::endl;

    std::cout << "Original C: ";
    displayMatrix(C);
    std::cout << std::endl;

    axpy(alpha, x, y);
    //daxpy(alpha, x, y);
    //dgemv(alpha, A, x, beta, y);
    gemv(alpha, A, x, beta, y);
    //dgemm(alpha, A, B, beta, C);
    gemm(alpha, A, B, beta, C);

    std::cout << "Final y: ";

    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << y[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Final C: ";
    displayMatrix(C);
    std::cout << std::endl;

    return 0;
}