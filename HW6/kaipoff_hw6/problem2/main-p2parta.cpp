#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <thread>
#include <chrono>
#include <stdexcept>
#include <mutex>
#include <mpi.h>

// Arc length integrand: sqrt(1 + (f'(x))^2)
double f(double x) {
    return std::sqrt(1 + std::pow((1 / x) - (x / 4.0), 2));
}

double sequential_riemann_sum(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = a + i * h;
        sum += f(x) * h;
    }
    return sum;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int ip, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &ip);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (argc != 2) {
        if (ip == 0)
            std::cerr << "Usage: mpirun -np <procs> ./program <npoints>\n";
        MPI_Finalize();
        return 1;
    }

    const int n = std::atoi(argv[1]);
    const double a = 1.0;
    const double b = 6.0;
    const double true_answer = 6.166759469;

    // Distribute work more evenly if n is not divisible by np
    int base_n = n / np;
    int remainder = n % np;
    int local_n = (ip < remainder) ? base_n + 1 : base_n;

    double h = (b - a) / n;
    double local_a = a + h * (ip * base_n + std::min(ip, remainder));
    double local_b = local_a + h * local_n;

    auto start = std::chrono::high_resolution_clock::now();
    double local_sum = sequential_riemann_sum(local_a, local_b, local_n);

    double global_sum;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    double elapsed_time = duration.count() * 1e-9;

    if (ip == 0) {
        double error = std::abs(global_sum - true_answer);
        std::cout << "Result: " << a << " to " << b << " = " << global_sum << std::endl;
        std::cout << "Parallel Sum:\t" << global_sum
                  << "\tTime: " << elapsed_time
                  << "\tError: " << error << std::endl;
    }

    MPI_Finalize();
    return 0;
}
