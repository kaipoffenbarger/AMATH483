#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <chrono>
#include <string>
#include <random>
#include "file_swaps.hpp"

int main(int argc, char* argv[]) {
    const int ntrial = 3;
    std::ofstream csv("swap_timing_results.csv");
    csv << "Dimension,RowSwapTime,ColSwapTime\n";

    std::mt19937 rng(std::random_device{}());

    for (int dim = 16; dim <= 8192; dim *= 2) {
        double rowTotalTime = 0.0;
        double colTotalTime = 0.0;
        std::string filename = "matrix_" + std::to_string(dim) + ".bin";

        for (int t = 0; t < ntrial; ++t) {
            int numRows = dim;
            int numCols = dim;

            // Generate the matrix
            std::vector<double> matrix(numRows * numCols);
            for (int col = 0; col < numCols; ++col)
                for (int row = 0; row < numRows; ++row)
                    matrix[col * numRows + row] = static_cast<double>(col * numRows + row);

            // Write the matrix to a file
            std::fstream file(filename, std::ios::out | std::ios::binary);
            file.write(reinterpret_cast<char*>(&matrix[0]), numRows * numCols * sizeof(double));
            file.close();

            // Get random indices i and j for row swapping
            std::uniform_int_distribution<int> dist(0, dim - 1);
            int i = dist(rng), j = dist(rng);
            while (i == j) j = dist(rng);

            // Row swap timing
            std::fstream fileToSwapRow(filename, std::ios::in | std::ios::out | std::ios::binary);
            auto startTimeRow = std::chrono::high_resolution_clock::now();
            swapRowsInFile(fileToSwapRow, numRows, numCols, i, j);
            auto endTimeRow = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> durationRow = endTimeRow - startTimeRow;
            rowTotalTime += durationRow.count();
            fileToSwapRow.close();

            // Column swap timing
            std::fstream fileToSwapCol(filename, std::ios::in | std::ios::out | std::ios::binary);
            auto startTimeCol = std::chrono::high_resolution_clock::now();
            swapColsInFile(fileToSwapCol, numRows, numCols, i, j);
            auto endTimeCol = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> durationCol = endTimeCol - startTimeCol;
            colTotalTime += durationCol.count();
            fileToSwapCol.close();

            // Delete the test file
            std::remove(filename.c_str());
        }

        double avgRowTime = rowTotalTime / ntrial;
        double avgColTime = colTotalTime / ntrial;
        csv << dim << "," << avgRowTime << "," << avgColTime << "\n";
        std::cout << "Dimension " << dim << " - Row: " << avgRowTime << "s, Col: " << avgColTime << "s\n";
    }

    csv.close();
    std::cout << "Results saved to swap_timing_results.csv\n";
    return 0;
}
