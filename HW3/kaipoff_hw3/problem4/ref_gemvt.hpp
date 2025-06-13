#include <iostream>
#include <cmath>
#include <vector>

#ifndef REF_GEMVT_HPP
#define REF_GEMVT_HPP

template <typename T>
void gemv(T a, const std::vector<std::vector<T>> &A, std::vector<T> &x, T b, std::vector<T> &y);

#endif