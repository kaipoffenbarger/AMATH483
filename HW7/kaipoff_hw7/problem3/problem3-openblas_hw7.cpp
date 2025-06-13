#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <cblas.h> // OpenBLAS header

void initialize_matrices(double* A, double* B, double* C, int n, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n * n; ++i) {
        A[i] = dist(gen);
        B[i] = dist(gen);
        C[i] = dist(gen);
    }
}

double measure_openblas_time(double* A, double* B, double* C, int n, int trials) {
    double total_time = 0.0;
    for (int t = 0; t < trials; ++t) {
        auto start = std::chrono::high_resolution_clock::now();
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 1.0, C, n);
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double>(end - start).count();
    }
    return total_time / trials;
}

int main() {
    // Initialize output file
    std::ofstream out("openblas_performance.csv");
    out << "n,OpenBLAS_GFLOPs\n";

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    const int trials = 3;

    // Test matrix sizes from 2 to 16384
    for (int k = 1; k <= 14; ++k) {
        int n = 1 << k; // n = 2^k
        std::cout << "Testing n = " << n << std::endl;

        // Allocate host matrices
        std::vector<double> A(n * n), B(n * n), C(n * n);
        initialize_matrices(A.data(), B.data(), C.data(), n, gen);

        // Measure OpenBLAS
        double openblas_time = measure_openblas_time(A.data(), B.data(), C.data(), n, trials);
        double openblas_gflops = (2.0 * n * n * n) / (openblas_time * 1e9);

        // Output results
        out << n << "," << openblas_gflops << "\n";
    }

    out.close();
    std::cout << "Results written to openblas_performance.csv\n";
    return 0;
}