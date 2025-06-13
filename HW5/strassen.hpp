// kaipoff@uw.edu
// AMATH 483-583
// strassen.hpp : Header file for Strassen implementation

#ifndef STRASSEN_HPP
#define STRASSEN_HPP

#include <vector>
#include <iostream>
#include <cmath>

// Function declarations
template <typename T>
std::vector<std::vector<T>> addMatrix(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);

template <typename T>
std::vector<std::vector<T>> subtractMatrix(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);

template <typename T>
std::vector<std::vector<T>> strassenMultiply(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);

template <typename T>
std::vector<std::vector<T>> classicMultiply(const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);

template <typename T>
std::vector<std::vector<T>> padMatrix(const std::vector<std::vector<T>> &A, int newSize);

template <typename T>
std::vector<std::vector<T>> unpadMatrix(const std::vector<std::vector<T>> &A, int originalSize);

int nextPowerOf2(int n);

template <typename T>
void printMatrix(const std::vector<std::vector<T>> &matrix);

#endif // STRASSEN_HPP