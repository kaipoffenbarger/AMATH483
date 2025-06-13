#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>

#include "strassen.hpp"

using Matrix = std::vector<std::vector<int>>;

int main() {

    Matrix A = {
        {1, 2, 4, 5, 6, 7},
        {3, 4, 6, 7, 8, 9},
        {1, 2, 4, 5, 6, 7},
        {3, 4, 6, 7, 8, 9},
        {1, 2, 4, 5, 6, 7},
        {3, 4, 6, 7, 8, 9}
    };
    Matrix B = {
        {1, 2, 4, 5, 6, 7},
        {3, 4, 6, 7, 8, 9},
        {1, 2, 4, 5, 6, 7},
        {3, 4, 6, 7, 8, 9},
        {1, 2, 4, 5, 6, 7},
        {3, 4, 6, 7, 8, 9}
    };

    Matrix result = strassenMultiply(A, B);

    std::cout << "Result of A x B:\n";
    printMatrix(result);


    return 0;
}