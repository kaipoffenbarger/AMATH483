#ifndef FILE_SWAPS_HPP
#define FILE_SWAPS_HPP

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

// Swap two rows in a matrix stored in column-major order in a binary file
void swapRowsInFile(std::fstream& file, int nRows, int nCols, int i, int j) {
    if (i == j) return;

    for (int col = 0; col < nCols; ++col) {
        std::streampos pos_i = static_cast<std::streampos>((col * nRows + i) * sizeof(double));
        std::streampos pos_j = static_cast<std::streampos>((col * nRows + j) * sizeof(double));

        double val_i, val_j;
        file.seekg(pos_i); file.read(reinterpret_cast<char*>(&val_i), sizeof(double));
        file.seekg(pos_j); file.read(reinterpret_cast<char*>(&val_j), sizeof(double));

        file.seekp(pos_i); file.write(reinterpret_cast<char*>(&val_j), sizeof(double));
        file.seekp(pos_j); file.write(reinterpret_cast<char*>(&val_i), sizeof(double));
    }
}

// Swap two columns in a matrix stored in column-major order in a binary file
void swapColsInFile(std::fstream& file, int nRows, int nCols, int i, int j) {
    if (i == j) return;

    for (int row = 0; row < nRows; ++row) {
        std::streampos pos_i = static_cast<std::streampos>((i * nRows + row) * sizeof(double));
        std::streampos pos_j = static_cast<std::streampos>((j * nRows + row) * sizeof(double));

        double val_i, val_j;
        file.seekg(pos_i); file.read(reinterpret_cast<char*>(&val_i), sizeof(double));
        file.seekg(pos_j); file.read(reinterpret_cast<char*>(&val_j), sizeof(double));

        file.seekp(pos_i); file.write(reinterpret_cast<char*>(&val_j), sizeof(double));
        file.seekp(pos_j); file.write(reinterpret_cast<char*>(&val_i), sizeof(double));
    }
}

#endif // FILE_SWAPS_HPP
