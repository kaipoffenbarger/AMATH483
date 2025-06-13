#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include "ref_daxpy.hpp"
#include "ref_dgemv.hpp"
#include "ref_dgemm.hpp"

int main()
{
    std::ofstream out("FLOPS_3.csv");
    out << "n,flops_per_sec\n";

    const int ntrials = 3;
    long double elapsed_time = 0.L;
    long double avg_time;

    std::mt19937 rng(5); // Random seed
    std::uniform_real_distribution<double> dist(1.0, 10.0); // Non-zero values

    for (int n = 2; n <= 512; ++n) {
        elapsed_time = 0.L; // Reset for this problem size

        for (int t = 0; t < ntrials; ++t) {
            double alpha = dist(rng);
            double beta = dist(rng);
            
            std::vector<std::vector<double>> A(n, std::vector<double>(n));
            std::vector<std::vector<double>> B(n, std::vector<double>(n));
            std::vector<std::vector<double>> C(n, std::vector<double>(n));

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    A[i][j] = dist(rng);
                    B[i][j] = dist(rng);
                    C[i][j] = dist(rng);
                }
            }

            auto start = std::chrono::high_resolution_clock::now();

            dgemm(alpha, A, B, beta, C);
            
            auto stop = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9); // in seconds

        }

        avg_time = elapsed_time / static_cast<long double>(ntrials);
        long double flops = static_cast<long double>((2*n*n*n) + (3*n*n));
        
        long double flops_per_sec = flops / avg_time;

        out << n << "," << flops_per_sec << "\n";
        std::cout << "Average Time = " << avg_time << ", n = " << n << ", Average FLOPs/sec = " << flops_per_sec << "\n";
    }

    out.close();
    return 0;
}
