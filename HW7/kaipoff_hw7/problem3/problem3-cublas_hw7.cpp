#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <cublas_v2.h>
#include <cuda_runtime.h>

void checkCudaError(cudaError_t err, const char* msg) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error: " << msg << ": " << cudaGetErrorString(err) << std::endl;
        exit(1);
    }
}

void initialize_matrices(double* A, double* B, double* C, int n, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n * n; ++i) {
        A[i] = dist(gen);
        B[i] = dist(gen);
        C[i] = dist(gen);
    }
}

double measure_cublas_time(cublasHandle_t handle, double* d_A, double* d_B, double* d_C, int n, int trials) {
    double total_time = 0.0;
    double alpha = 1.0, beta = 1.0;
    for (int t = 0; t < trials; ++t) {
        auto start = std::chrono::high_resolution_clock::now();
        cublasStatus_t status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &alpha, d_A, n, d_B, n, &beta, d_C, n);
        if (status != CUBLAS_STATUS_SUCCESS) {
            std::cerr << "cuBLAS Error: DGEMM failed" << std::endl;
            exit(1);
        }
        checkCudaError(cudaDeviceSynchronize(), "Device synchronization failed");
        auto end = std::chrono::high_resolution_clock::now();
        total_time += std::chrono::duration<double>(end - start).count();
    }
    return total_time / trials;
}

int main() {
    // Initialize output file
    std::ofstream out("cublas_performance.csv");
    out << "n,cuBLAS_GFLOPs\n";

    // Initialize cuBLAS
    cublasHandle_t handle;
    cublasStatus_t status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        std::cerr << "cuBLAS Error: Initialization failed" << std::endl;
        return 1;
    }

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

        // Allocate device matrices
        double *d_A, *d_B, *d_C;
        checkCudaError(cudaMalloc(&d_A, n * n * sizeof(double)), "cudaMalloc d_A");
        checkCudaError(cudaMalloc(&d_B, n * n * sizeof(double)), "cudaMalloc d_B");
        checkCudaError(cudaMalloc(&d_C, n * n * sizeof(double)), "cudaMalloc d_C");
        checkCudaError(cudaMemcpy(d_A, A.data(), n * n * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy d_A");
        checkCudaError(cudaMemcpy(d_B, B.data(), n * n * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy d_B");
        checkCudaError(cudaMemcpy(d_C, C.data(), n * n * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy d_C");

        // Measure cuBLAS
        double cublas_time = measure_cublas_time(handle, d_A, d_B, d_C, n, trials);
        double cublas_gflops = (2.0 * n * n * n) / (cublas_time * 1e9);

        // Output results
        out << n << "," << cublas_gflops << "\n";

        // Free device memory
        checkCudaError(cudaFree(d_A), "cudaFree d_A");
        checkCudaError(cudaFree(d_B), "cudaFree d_B");
        checkCudaError(cudaFree(d_C), "cudaFree d_C");
    }

    // Cleanup
    cublasDestroy(handle);
    out.close();
    std::cout << "Results written to cublas_performance.csv\n";
    return 0;
}