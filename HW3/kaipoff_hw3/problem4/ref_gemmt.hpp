#include <iostream>
#include <cmath>
#include <vector>

#ifndef REF_GEMMT_HPP
#define REF_GEMMT_HPP

template <typename T>
void gemm(T a, const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B, T b, std::vector<std::vector<T>> &C);

#endif