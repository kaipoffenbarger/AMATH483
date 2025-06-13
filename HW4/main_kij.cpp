#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include "kij_loop.hpp"

int main()
{   

    // DOUBLE
    std::cout << "Calculating Doubles [KIJ]" << std::endl;

    std::ofstream out("kij_double_O3.csv");
    out << "n,flops_per_sec\n";

    const int ntrials = 3;
    long double elapsed_time = 0.L;
    long double avg_time;

    std::mt19937 rng(5); // Random seed
    std::uniform_real_distribution<double> dis(1.0, 10.0); // Non-zero values

    for (int n = 2; n <= 512; ++n) {
        elapsed_time = 0.L; // Reset for this problem size
        int total_vals = n*n;
        double alpha = dis(rng);
        double beta = dis(rng);
        std::vector<double> A(total_vals), B(total_vals), C(total_vals);
        
        for (int i = 0; i < total_vals; ++i) {
            A[i] = dis(rng);
            B[i] = dis(rng);
            C[i] = dis(rng);
        }

        for (int t = 0; t < ntrials; ++t) {

            auto start = std::chrono::high_resolution_clock::now();

            mm_kij(alpha, A, B, beta, C, n, n, n);
            
            auto stop = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9); // in seconds

        }

        avg_time = elapsed_time / static_cast<long double>(ntrials);
        long double flops = static_cast<long double>((3*n*n*n) + (n*n));
        long double flops_per_sec = flops / avg_time;

        out << n << "," << flops_per_sec << "\n";
        std::cout << "[DOUBLE, KIJ, O3] Average Time = " << avg_time << ", n = " << n << ", Average FLOPs/sec = " << flops_per_sec << "\n";
    }

    out.close();

    // FLOAT
    std::cout << "Now onto floats [KIJ]" << std::endl;

    std::ofstream out2("kij_float_O3.csv");
    out2 << "n,flops_per_sec\n";

    elapsed_time = 0.L;

    std::mt19937 rng2(6); // Random seed
    std::uniform_real_distribution<float> dist(1.0, 10.0); // Non-zero values

    for (int n = 2; n <= 512; ++n) {
        elapsed_time = 0.L; // Reset for this problem size
        int total_vals = n*n;
        float alpha = dist(rng2);
        float beta = dist(rng2);
        std::vector<float> A(total_vals), B(total_vals), C(total_vals);
        
        for (int i = 0; i < total_vals; ++i) {
            A[i] = dist(rng2);
            B[i] = dist(rng2);
            C[i] = dist(rng2);
        }

        for (int t = 0; t < ntrials; ++t) {

            auto start = std::chrono::high_resolution_clock::now();

            mm_kij(alpha, A, B, beta, C, n, n, n);
            
            auto stop = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time += (duration.count() * 1.e-9); // in seconds

        }

        avg_time = elapsed_time / static_cast<long double>(ntrials);
        long double flops = static_cast<long double>((3*n*n*n) + (n*n));
        long double flops_per_sec = flops / avg_time;

        out2 << n << "," << flops_per_sec << "\n";
        std::cout << "[FLOAT, KIJ, O3] Average Time = " << avg_time << ", n = " << n << ", Average FLOPs/sec = " << flops_per_sec << "\n";
    }

    out2.close();

    return 0;
}