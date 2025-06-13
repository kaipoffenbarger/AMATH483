#include <iostream>
#include <cmath>
#include <vector>

#ifndef REF_DGEMV_HPP
#define REF_DGEMV_HPP

void dgemv(double a, const std::vector<std::vector<double>> &A, std::vector<double> &x, double b, std::vector<double> &y);

#endif