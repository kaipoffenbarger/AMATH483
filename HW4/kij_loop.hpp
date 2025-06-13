#ifndef KIJ_LOOP_HPP
#define KIJ_LOOP_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

template <typename T>
void mm_kij (T a, const std::vector<T> &A, const std::vector<T> &B, T b, std::vector<T> &C, int m, int p, int n);

#endif // KIJ_LOOP_HPP