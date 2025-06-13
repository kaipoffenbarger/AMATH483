#include <fftw3.h>
#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <cmath>
#include <fstream>

int main() {
    std::ofstream out("fftw_perf.csv");
    out << "n,flops\n";

    const int ntrial = 3;

    for (int n = 16; n <= 256; n *= 2) {
        size_t N = n * n * n;
        double total_time = 0.0;

        for (int trial = 0; trial < ntrial; ++trial) {
            fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

            for (size_t i = 0; i < N; ++i)
                in[i][0] = cos(2*M_PI*i/n), in[i][1] = sin(2*M_PI*i/n);

            auto plan_fwd = fftw_plan_dft_3d(n, n, n, in, out, FFTW_FORWARD, FFTW_MEASURE);
            auto plan_bwd = fftw_plan_dft_3d(n, n, n, out, in, FFTW_BACKWARD, FFTW_MEASURE);

            auto start = std::chrono::high_resolution_clock::now();
            fftw_execute(plan_fwd);
            // simulate multiply by ik in Fourier space
            fftw_execute(plan_bwd);
            auto end = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> elapsed = end - start;
            total_time += elapsed.count();

            fftw_destroy_plan(plan_fwd);
            fftw_destroy_plan(plan_bwd);
            fftw_free(in); fftw_free(out);
        }

        double avg_time = total_time / ntrial;
        double flops = (30 * std::pow(n,3) * std::log2(n) + 9 * std::pow(n,3)) / avg_time;
        out << n << "," << flops << "\n";
    }

    return 0;
}
