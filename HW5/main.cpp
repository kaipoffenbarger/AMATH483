#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>

#include "strassen.hpp"
using Matrix = std::vector<std::vector<double>>;

long long strassenFLOPs(int n, int base_case_size = 32) {
    if (n <= base_case_size) {
        return 2LL * n * n * n; // Classic triple-loop
    }

    int half = n / 2;

    long long recursive = 7LL * strassenFLOPs(half, base_case_size);
    long long adds = 18LL * half * half;

    return recursive + adds;
}

int main() {
    
    std::ofstream out("strassen.csv");
    out << "n,flops_per_sec\n";

    const int ntrials = 3;
    long double elapsed_time = 0.L;
    long double avg_time;

    std::mt19937 rng(5); // Random seed
    std::uniform_real_distribution<double> dist(1.0, 10.0); // Non-zero values

    for (int n = 2; n <= 512; n += 2) {
        elapsed_time = 0.L; // Reset for this problem size

        Matrix A(n, std::vector<double>(n));
        Matrix B(n, std::vector<double>(n));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = dist(rng);
                B[i][j] = dist(rng);
            }
        }

        for (int t = 0; t < ntrials; ++t) {

            auto start = std::chrono::steady_clock::now();

            strassenMultiply(A, B);
            
            auto stop = std::chrono::steady_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9); // in seconds

        }

        if (elapsed_time < 0) {
            std::cout << "Elapsed Time = " << elapsed_time << ", n = " << n << ". Average FLOPs/sec not recorded because its negative." << "\n";
            continue;
        }
        avg_time = elapsed_time / static_cast<long double>(ntrials);
        long double flops = static_cast<long double>(strassenFLOPs(n));
        long double flops_per_sec = flops / avg_time;

        out << n << "," << flops_per_sec << "\n";
        std::cout << "Average Time = " << avg_time << ", n = " << n << ", Average FLOPs/sec = " << flops_per_sec << "\n";
    }

    out.close();
    return 0;

}