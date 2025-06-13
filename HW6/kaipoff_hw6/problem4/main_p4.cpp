#include <mpi.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include "my_broadcast.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            std::cerr << "Usage: ./xmpi_bandwidth <custom|mpi> <repeats>\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string mode = argv[1];
    int repeats = std::atoi(argv[2]);

    std::vector<size_t> sizes;
    for (size_t s = 8; s <= 268435456; s *= 2) // 8B to 256MB
        sizes.push_back(s);

    if (rank == 0) {
        std::cout << "MessageSize(Bytes),Bandwidth(Bytes/sec)\n";
    }

    for (size_t size_bytes : sizes) {
        std::vector<char> buffer(size_bytes, rank);  // Initialize buffer with rank for validation
        MPI_Barrier(MPI_COMM_WORLD);

        double total_time = 0.0;
        for (int i = 0; i < repeats; ++i) {
            MPI_Barrier(MPI_COMM_WORLD);
            auto start = std::chrono::high_resolution_clock::now();

            if (mode == "custom") {
                my_broadcast(buffer.data(), size_bytes, 0, MPI_COMM_WORLD);
            } else {
                MPI_Bcast(buffer.data(), size_bytes, MPI_BYTE, 0, MPI_COMM_WORLD);
            }

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            total_time += elapsed.count();
        }

        double avg_time = total_time / repeats;
        double bandwidth = size_bytes / avg_time;

        if (rank == 0) {
            std::cout << size_bytes << "," << std::fixed << std::setprecision(2) << bandwidth << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}