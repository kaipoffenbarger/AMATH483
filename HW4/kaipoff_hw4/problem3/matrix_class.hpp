// AMATH 483-583 row major Matrix class template starter
// write the methods for:
// transpose
// infinityNorm
// operator*
// operator+

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>

template <typename T>
class Matrix
{
public:
    Matrix(int numRows, int numCols)
        : num_rows(numRows), num_cols(numCols), data(numRows * numCols) {}

    T &operator()(int i, int j)
    {
        return data[i * num_cols + j];
    }

    const T &operator()(int i, int j) const
    {
        return data[i * num_cols + j];
    }

    Matrix<T> operator*(const Matrix<T> &other) const
    {
        // Dimension check
        if (num_cols != other.numRows()) {
            throw std::invalid_argument("[Multiplication] Column dimension doesn't match row dimension of input matrix!");
        }

        Matrix<T> result = Matrix(num_rows, other.numCols());

        for (int i = 0; i < num_rows; ++i) {
            for (int j = 0; j < other.numCols(); ++j) {
                for (int k = 0; k < other.numRows(); ++k) {
                    result(i, j) += (*this)(i, k) * other(k, j); //value from B is fixed in this loop
                }
            }
        }

        return result;

    }

    Matrix<T> operator+(const Matrix<T> &other) const;

    Matrix<T> transpose() const
    {
        Matrix<T> temp = Matrix(num_cols, num_rows);

        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                temp(j, i) = (*this)(i, j);
            }
        }
        
        return temp;
    }

    int numRows() const
    {
        return num_rows;
    }

    int numCols() const
    {
        return num_cols;
    }

    T infinityNorm() const
    {
        T norm = 0;
        std::vector<T> rowsums(num_rows);
        for (int i = 0; i < num_rows; i++) {
            rowsums[i] = 0;
            for (int j = 0; j < num_cols; j++){
                rowsums[i] += std::abs((*this)(i, j));
            }
        }
        norm = *std::max_element(rowsums.begin(), rowsums.end());
        return norm;
    }

private:
    int num_rows;
    int num_cols;
    std::vector<T> data;
};

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const
{
    // check for errors in Matrix dimensions
    if (num_cols != other.numCols()) {
        throw std::invalid_argument("[Addition] Column dimension doesn't match row dimension of input matrix!");
    }
    if (num_rows != other.numRows()) {
        throw std::invalid_argument("[Addition] Row dimension doesn't match row dimension of input matrix!");
    }

    Matrix<T> result(num_rows, num_cols);

    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            result(i, j) = (*this)(i, j) + other(i, j);
        }
    }

    return result;
}
