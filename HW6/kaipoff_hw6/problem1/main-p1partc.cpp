#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <thread>
#include <chrono>
#include <stdexcept>
#include <mutex>

double f(double x) {
    return std::sqrt(1 + std::pow((1 / x) - (x / 4), 2));
}

double riemann_sum(double a, double b, int n) {
    double h = (b-a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x = a + i*h;
        sum += f(x);
    }
    sum *= h;
    return sum;
}

void compute_partial_sum(double &sum, std::mutex &sum_mutex, double a, double h, int n, int num_threads, int i) {
    double partial_sum = 0.0;
    int start = (n / num_threads) * i;
    int end = (i == num_threads - 1) ? n: (n/num_threads) * (i + 1);

    for (int j = start; j < end; ++j) {
        double x = a + j * h;
        partial_sum += f(x);
    }
    partial_sum *= h;

    //lock mutex then update sum
    sum_mutex.lock();
    sum += partial_sum;
    sum_mutex.unlock();
}

double parrallel_riemann_sum(double a, double b, int n, int num_threads) {
    double h = (b-a) /n;
    double sum = 0.0;
    std::vector<std::thread> threads(num_threads);

    std::mutex sum_mutex;

    //compute sums (thread spawning)
    for (int i=0; i < num_threads; ++i) {
        threads[i] = std::thread(compute_partial_sum, std::ref(sum), std::ref(sum_mutex), a, h, n, num_threads, i);
    }

    //wait for threads to finish tasks
    for (int i =0; i < num_threads; ++i) {
        threads[i].join();
    }

    return sum;
}


int main(int argc, char *argv[]) {

    const double a = 1.0;
    const double b = 6.0;
    const int n = 100000000; // change for part d
    long double elapsed_time = 0.L;

    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " nthreads " << "npoints" << std::endl;
        return 1;
    }

    std::vector<int> thread_counts = {1, 2, 4, 8, 16}; //was originally trying it with a for loop
    const int num_threads = std::atoi(argv[1]);
    const int npoints = std::atoi(argv[2]);

    //Part b and c (python code for plotting)

    elapsed_time = 0.L;
    start = std::chrono::high_resolution_clock::now();
    double par_sum = parrallel_riemann_sum(a, b, npoints, num_threads);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    elapsed_time = (duration.count() * 1.e-9);
    std::cout << "Parallel:\t" << par_sum << "\t" << "time: " << elapsed_time << std::endl;

    return 0;
}