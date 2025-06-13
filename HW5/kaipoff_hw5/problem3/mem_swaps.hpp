#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

void swapRows(std::vector<double> &matrix, int nRows, int nCols, int i, int j){

    //check if matrix is empty
    if (matrix.empty()) {
        throw std::invalid_argument("Matrix is empty!");
    }

    //make sure there are enough rows
    if (nRows <= 1) {
        throw std::invalid_argument("Matrix needs at least two rows for row swap!");
    }

    //using zero indexing for i and j
    if (i >= nRows || i < 0) {
        throw std::invalid_argument("Row index i out of bounds.");
    }

    if (j >= nRows || j < 0) {
        throw std::invalid_argument("Row index j out of bounds.");
    }

    //swap rows
    for (int col = 0; col < nCols; ++col) {
        std::swap(matrix[i + col * nRows], matrix[j + col * nRows]);
    }


}


void swapCols(std::vector<double> &matrix, int nRows, int nCols, int i, int j){
    //check if matrix is empty
    if (matrix.empty()) {
        throw std::invalid_argument("Matrix is empty!");
    }

    //make sure there are enough rows
    if (nCols <= 1) {
        throw std::invalid_argument("Matrix needs at least two rows for row swap!");
    }

    //using zero indexing for i and j
    if (i >= nCols || i < 0) {
        throw std::invalid_argument("Column index i out of bounds.");
    }

    if (j >= nCols || j < 0) {
        throw std::invalid_argument("Column index j out of bounds.");
    }

    //swap cols
    for (int row = 0; row < nRows; ++row) {
        std::swap(matrix[row + i * nRows], matrix[row + j * nRows]);
    }
}