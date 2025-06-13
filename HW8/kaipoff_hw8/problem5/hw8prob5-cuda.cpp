#include <cufft.h>
#include <cuda_runtime.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        std::cerr << "CUDA Error: " << cudaGetErrorString(code) << " at " << file << ":" << line << std::endl;
        exit(code);
    }
}

__global__ void multiply_by_ik(cufftDoubleComplex* data, int n, int dir, double L) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int N = n * n * n;
    if (idx >= N) return;

    int i = idx % n;
    int j = (idx / n) % n;
    int k = (idx / (n * n));

    int kx = (i <= n/2) ? i : i - n;
    int ky = (j <= n/2) ? j : j - n;
    int kz = (k <= n/2) ? k : k - n;

    int kdir = (dir == 0) ? kx : (dir == 1) ? ky : kz;
    double phase = 2 * M_PI * kdir / L;

    double real = data[idx].x;
    double imag = data[idx].y;

    data[idx].x = -phase * imag;
    data[idx].y =  phase * real;
}

int main() {
    std::ofstream out("cufft_perf.csv");
    out << "n,flops\n";

    const int ntrial = 3;
    const int stride = 2;

    for (int n = 16; n <= 256; n *= stride) {
        int N = n * n * n;
        size_t size = sizeof(cufftDoubleComplex) * N;

        double total_time = 0.0;

        for (int trial = 0; trial < ntrial; ++trial) {
            cufftDoubleComplex* h_data = new cufftDoubleComplex[N];
            for (int i = 0; i < N; ++i) {
                h_data[i].x = cos(2 * M_PI * i / n);
                h_data[i].y = sin(2 * M_PI * i / n);
            }

            cufftDoubleComplex *d_data;
            CUDA_CHECK(cudaMalloc(&d_data, size));
            CUDA_CHECK(cudaMemcpy(d_data, h_data, size, cudaMemcpyHostToDevice));

            cufftHandle plan;
            cufftPlan3d(&plan, n, n, n, CUFFT_Z2Z);

            cudaEvent_t start, stop;
            cudaEventCreate(&start);
            cudaEventCreate(&stop);

            cudaEventRecord(start);

            for (int dir = 0; dir < 3; ++dir) {
                cufftExecZ2Z(plan, d_data, d_data, CUFFT_FORWARD);
                multiply_by_ik<<<(N + 255)/256, 256>>>(d_data, n, dir, n);
                cufftExecZ2Z(plan, d_data, d_data, CUFFT_INVERSE);
            }

            cudaEventRecord(stop);
            cudaEventSynchronize(stop);

            float milliseconds = 0;
            cudaEventElapsedTime(&milliseconds, start, stop);

            total_time += milliseconds / 1000.0;

            delete[] h_data;
            cufftDestroy(plan);
            cudaFree(d_data);
            cudaEventDestroy(start);
            cudaEventDestroy(stop);
        }

        double avg_time = total_time / ntrial;
        double flops = (30.0 * N * std::log2(n) + 9.0 * N) / avg_time;
        out << n << "," << flops << "\n";
    }

    return 0;
}
