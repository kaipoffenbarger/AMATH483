#include <iostream>
#include <vector>
#include <chrono>
#include <cblas.h>
#include <fstream>
#include <iomanip>

const int N_START = 2;
const int N_END = 4096;
const int NTRIAL = 3;
const int STRIDE = 2;

void benchmark_daxpy(std::ostream& out) {
    out << "# L1 daxpy: n MFLOPs\n";
    for (int n = N_START; n <= N_END; n *= STRIDE) {
        std::vector<double> x(n, 1.0), y(n, 2.0);
        double alpha = 2.5;
        double total_time = 0.0;

        for (int trial = 0; trial < NTRIAL; ++trial) {
            std::fill(y.begin(), y.end(), 2.0);
            auto start = std::chrono::high_resolution_clock::now();
            cblas_daxpy(n, alpha, x.data(), 1, y.data(), 1);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double>(end - start).count();
        }

        double avg_time = total_time / NTRIAL;
        double flops = 2.0 * n;
        double mflops = (flops / avg_time) / 1e6;
        out << n << " " << std::fixed << std::setprecision(2) << mflops << "\n";
    }
    out << "\n";
}

void benchmark_dgemv(std::ostream& out) {
    out << "# L2 dgemv: n MFLOPs\n";
    for (int n = N_START; n <= N_END; n *= STRIDE) {
        std::vector<double> A(n * n, 1.0), x(n, 1.0), y(n, 0.0);
        double alpha = 1.0, beta = 0.0;
        double total_time = 0.0;

        for (int trial = 0; trial < NTRIAL; ++trial) {
            std::fill(y.begin(), y.end(), 0.0);
            auto start = std::chrono::high_resolution_clock::now();
            cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, A.data(), n,
                        x.data(), 1, beta, y.data(), 1);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double>(end - start).count();
        }

        double avg_time = total_time / NTRIAL;
        double flops = 2.0 * n * n;
        double mflops = (flops / avg_time) / 1e6;
        out << n << " " << std::fixed << std::setprecision(2) << mflops << "\n";
    }
    out << "\n";
}

void benchmark_dgemm(std::ostream& out) {
    out << "# L3 dgemm: n MFLOPs\n";
    for (int n = N_START; n <= N_END; n *= STRIDE) {
        std::vector<double> A(n * n, 1.0), B(n * n, 1.0), C(n * n, 0.0);
        double alpha = 1.0, beta = 0.0;
        double total_time = 0.0;

        for (int trial = 0; trial < NTRIAL; ++trial) {
            std::fill(C.begin(), C.end(), 0.0);
            auto start = std::chrono::high_resolution_clock::now();
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        n, n, n, alpha, A.data(), n, B.data(), n, beta, C.data(), n);
            auto end = std::chrono::high_resolution_clock::now();
            total_time += std::chrono::duration<double>(end - start).count();
        }

        double avg_time = total_time / NTRIAL;
        double flops = 2.0 * n * n * n;
        double mflops = (flops / avg_time) / 1e6;
        out << n << " " << std::fixed << std::setprecision(2) << mflops << "\n";
    }
    out << "\n";
}

int main() {
    std::ofstream outfile("openblas_perf.txt");

    if (!outfile.is_open()) {
        std::cerr << "Error opening output file!\n";
        return 1;
    }

    benchmark_daxpy(outfile);
    benchmark_dgemv(outfile);
    benchmark_dgemm(outfile);

    outfile.close();
    std::cout << "Benchmarking complete. Results saved to openblas_perf.txt\n";
    return 0;
}
