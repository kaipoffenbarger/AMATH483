#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <utility>

#include "mem_swaps.hpp"

std::pair<int, int> getRandomIndices(int n) {
    int i = std::rand() % n;
    int j = std::rand() % ( n - 1);
    if (j >= i) {
    j++;
    }
    return std::make_pair(i, j);
}

void print_matrix(const std::vector<double>& data, size_t rows, size_t cols) {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << data[i + j * rows] << " ";
        }
        std::cout << "\n";
    }
}

int main() {

    const int ntrials = 3;
    std::ofstream out("row_swap.csv");
    out << "n,average_time\n";
    std::ofstream out2("col_swap.csv");
    out2 << "n,average_time\n";

    long double elapsed_time_row = 0.L;
    long double avg_time_row;

    long double elapsed_time_col = 0.L;
    long double avg_time_col;

    std::mt19937 rng(5); // Random seed
    std::uniform_real_distribution<double> dist(1.0, 10.0); // Non-zero values

    //iterate over each n
    for (int n = 16; n <= 4096; n *= 2) {

        elapsed_time_row = 0.L; // Reset for this problem size
        elapsed_time_col = 0.L;

        //initalize matrix
        std::vector<double> matrix(n*n);
        for (int i = 0; i < (n*n); i++) {
            matrix[i] = dist(rng);
        }

        //execute trials
        for (int i = 0; i < ntrials; i++) {

            //ROW: get indicies
            std::pair<int, int> rowIndices = getRandomIndices(n);
            int j = rowIndices.first;
            int k = rowIndices.second;

            auto start = std::chrono::steady_clock::now();

            swapRows(matrix, n, n, j, k);
            
            auto stop = std::chrono::steady_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
            elapsed_time_row += (duration.count() * 1.e-9); // in seconds



            //COL: get indicies
            std::pair<int, int> colIndices = getRandomIndices(n);
            int l = colIndices.first;
            int m = colIndices.second;

            auto start2 = std::chrono::steady_clock::now();

            swapCols(matrix, n, n, l, m);
            
            auto stop2 = std::chrono::steady_clock::now();

            auto duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(stop2 - start2);
            elapsed_time_col += (duration2.count() * 1.e-9); // in seconds

        }

        if (elapsed_time_row < 0) {
            std::cout << "Elapsed Time Row = " << elapsed_time_row << ", n = " << n << ". Time not recorded because its negative." << "\n";
            continue;
        }

        if (elapsed_time_col < 0) {
            std::cout << "Elapsed Time Col = " << elapsed_time_row << ", n = " << n << ". Time not recorded because its negative." << "\n";
            continue;
        }

        avg_time_row = elapsed_time_row / static_cast<long double>(ntrials);
        avg_time_col = elapsed_time_col / static_cast<long double>(ntrials);

        out << n << "," << avg_time_row << "\n";
        std::cout << "[ROW] Average Time = " << avg_time_row << ", n = " << n << "\n";

        out2 << n << "," << avg_time_col << "\n";
        std::cout << "[COL] Average Time = " << avg_time_col << ", n = " << n << "\n";

    }

    out.close();
    out2.close();
    return 0;
}