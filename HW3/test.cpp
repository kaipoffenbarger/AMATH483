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
    std::ofstream out("FLOPS_1.csv");
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
            std::vector<double> x(n), y(n), y_original(n);

            for (int i = 0; i < n; ++i) {
                x[i] = dist(rng);
                y[i] = dist(rng);
                y_original[i] = y[i];
            }

            auto start = std::chrono::high_resolution_clock::now();

            daxpy(alpha, x, y);
            
            auto stop = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9); // in seconds

        }

        avg_time = elapsed_time / static_cast<long double>(ntrials);
        long double flops = static_cast<long double>(2 * n);
        long double flops_per_sec = flops / avg_time;

        out << n << "," << flops_per_sec << "\n";
        std::cout << "Average Time = " << avg_time << ", n = " << n << ", Average FLOPs/sec = " << flops_per_sec << "\n";
    }

    out.close();
    return 0;
}
