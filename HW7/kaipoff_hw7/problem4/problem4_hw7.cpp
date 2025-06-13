#include <iostream>
#include <fstream>
#include <vector>
#include <cuda_runtime.h>
#include <cuda.h>
#include <iomanip>

void checkCudaError(cudaError_t err, const char* msg) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error: " << msg << ": " << cudaGetErrorString(err) << std::endl;
        exit(1);
    }
}

double measure_bandwidth(size_t size, int trials, cudaMemcpyKind kind) {
    double* h_data;
    double* d_data;
    cudaEvent_t start, stop;
    checkCudaError(cudaMallocHost(&h_data, size), "Host allocation");
    checkCudaError(cudaMalloc(&d_data, size), "Device allocation");
    checkCudaError(cudaEventCreate(&start), "Event create start");
    checkCudaError(cudaEventCreate(&stop), "Event create stop");

    // Initialize host data
    for (size_t i = 0; i < size / sizeof(double); ++i) {
        h_data[i] = static_cast<double>(i);
    }

    double total_time = 0.0;
    for (int t = 0; t < trials; ++t) {
        checkCudaError(cudaEventRecord(start, 0), "Event record start");
        if (kind == cudaMemcpyHostToDevice) {
            checkCudaError(cudaMemcpy(d_data, h_data, size, cudaMemcpyHostToDevice), "H2D memcpy");
        } else {
            checkCudaError(cudaMemcpy(h_data, d_data, size, cudaMemcpyDeviceToHost), "D2H memcpy");
        }
        checkCudaError(cudaEventRecord(stop, 0), "Event record stop");
        checkCudaError(cudaEventSynchronize(stop), "Event synchronize");
        float ms;
        checkCudaError(cudaEventElapsedTime(&ms, start, stop), "Event elapsed time");
        total_time += ms / 1000.0; // Convert to seconds
    }

    double avg_time = total_time / trials;
    double bandwidth = (size / (avg_time * 1e9)); // GB/s
    checkCudaError(cudaFreeHost(h_data), "Free host");
    checkCudaError(cudaFree(d_data), "Free device");
    checkCudaError(cudaEventDestroy(start), "Destroy start event");
    checkCudaError(cudaEventDestroy(stop), "Destroy stop event");
    return bandwidth;
}

int main() {
    std::ofstream out("bandwidth.csv");
    out << "BufferSize_Bytes,H2D_Bandwidth_GBs,D2H_Bandwidth_GBs\n";

    const int trials = 10;
    for (int k = 0; k <= 31; ++k) {
        size_t size = 1ULL << k; // 2^k bytes
        std::cout << "Testing buffer size: " << size << " bytes\n";

        double h2d_bandwidth = measure_bandwidth(size, trials, cudaMemcpyHostToDevice);
        double d2h_bandwidth = measure_bandwidth(size, trials, cudaMemcpyDeviceToHost);

        out << size << "," << std::fixed << std::setprecision(3) << h2d_bandwidth << ","
            << d2h_bandwidth << "\n";
    }

    out.close();
    std::cout << "Results written to bandwidth.csv\n";
    return 0;
}