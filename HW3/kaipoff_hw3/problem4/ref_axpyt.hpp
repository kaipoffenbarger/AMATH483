#include <iostream>
#include <cmath>
#include <vector>

#ifndef REF_AXPYT_HPP
#define REF_AXPYT_HPP

template<typename T>
void axpy(T a, const std::vector<T> &x, std::vector<T> &y);

#endif