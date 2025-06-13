#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

template <typename T>
void mm_kij(T a, const std::vector<T>& A, const std::vector<T>& B, T b, std::vector<T>& C, int m, int p, int n) {
    // Scale C by beta first
    for (int i = 0; i < m * n; ++i) {
        C[i] *= b;
    }

    // Perform matrix multiplication: C += alpha * A * B
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < m; ++j) {
                C[i * n + j] += a * A[i * p + k] * B[k * n + j]; //value from B is fixed in this loop
            }
        }
    }
}

template void mm_kij(float a, const std::vector<float> &A, const std::vector<float> &B, float b, std::vector<float> &C, int, int, int);
template void mm_kij(double a, const std::vector<double> &A, const std::vector<double> &B, double b, std::vector<double> &C, int, int, int);