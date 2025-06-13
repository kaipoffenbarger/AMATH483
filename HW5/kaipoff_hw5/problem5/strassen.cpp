// kaipoff@uw.edu
// AMATH 483-583
// strassen.cpp : starter code for Strassen implementation

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template <typename T>
vector<vector<T>> addMatrix(const vector<vector<T>> &A, const vector<vector<T>> &B)
{
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> C(n, vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
    return C;
}

template <typename T>
vector<vector<T>> subtractMatrix(const vector<vector<T>> &A, const vector<vector<T>> &B)
{
    int n = A.size();
    int m = A[0].size();
    vector<vector<T>> C(n, vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

// Helper: Get next power of 2
int nextPowerOf2(int n) {
    if (n == 0) return 1;
    return std::pow(2, std::ceil(std::log2(n)));
}

// Helper: Pad matrix with zeros to size newSize x newSize
template <typename T>
std::vector<std::vector<T>> padMatrix(const std::vector<std::vector<T>> &A, int newSize) {
    int oldRows = A.size();
    int oldCols = A[0].size();
    std::vector<std::vector<T>> padded(newSize, std::vector<T>(newSize, 0));
    for (int i = 0; i < oldRows; ++i)
        for (int j = 0; j < oldCols; ++j)
            padded[i][j] = A[i][j];
    return padded;
}

// Helper: Unpad matrix back to originalSize x originalSize
template <typename T>
std::vector<std::vector<T>> unpadMatrix(const std::vector<std::vector<T>> &A, int originalSize) {
    std::vector<std::vector<T>> result(originalSize, std::vector<T>(originalSize));
    for (int i = 0; i < originalSize; ++i)
        for (int j = 0; j < originalSize; ++j)
            result[i][j] = A[i][j];
    return result;
}

template <typename T>
vector<vector<T>> classicMultiply(const vector<vector<T>> &A, const vector<vector<T>> &B) {
    int n = A.size();
    vector<vector<T>> C(n, vector<T>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

template <typename T>
vector<vector<T>> strassenMultiply(const vector<vector<T>> &A, const vector<vector<T>> &B)
{
    int n = A.size();
    const int threshold = 32;

    bool padded = false;
    int newSize = n;

    if ((n & (n - 1)) != 0) { // Not a power of two
        newSize = std::pow(2, std::ceil(std::log2(n)));
        padded = true;
    }

    vector<vector<T>> A_pad = padded ? padMatrix(A, newSize) : A;
    vector<vector<T>> B_pad = padded ? padMatrix(B, newSize) : B;

    if (newSize <= threshold) {
        auto C_pad = classicMultiply(A_pad, B_pad);
        return padded ? unpadMatrix(C_pad, n) : C_pad;
    }

    int new_dim = newSize / 2;

    vector<vector<T>> A11(new_dim, vector<T>(new_dim));
    vector<vector<T>> A12(new_dim, vector<T>(new_dim));
    vector<vector<T>> A21(new_dim, vector<T>(new_dim));
    vector<vector<T>> A22(new_dim, vector<T>(new_dim));
    vector<vector<T>> B11(new_dim, vector<T>(new_dim));
    vector<vector<T>> B12(new_dim, vector<T>(new_dim));
    vector<vector<T>> B21(new_dim, vector<T>(new_dim));
    vector<vector<T>> B22(new_dim, vector<T>(new_dim));

    for (int i = 0; i < new_dim; ++i) {
        for (int j = 0; j < new_dim; ++j) {
            A11[i][j] = A_pad[i][j];
            A12[i][j] = A_pad[i][j + new_dim];
            A21[i][j] = A_pad[i + new_dim][j];
            A22[i][j] = A_pad[i + new_dim][j + new_dim];
            B11[i][j] = B_pad[i][j];
            B12[i][j] = B_pad[i][j + new_dim];
            B21[i][j] = B_pad[i + new_dim][j];
            B22[i][j] = B_pad[i + new_dim][j + new_dim];
        }
    }

    auto P1 = strassenMultiply(addMatrix(A11, A22), addMatrix(B11, B22));
    auto P2 = strassenMultiply(addMatrix(A21, A22), B11);
    auto P3 = strassenMultiply(A11, subtractMatrix(B12, B22));
    auto P4 = strassenMultiply(A22, subtractMatrix(B21, B11));
    auto P5 = strassenMultiply(addMatrix(A11, A12), B22);
    auto P6 = strassenMultiply(subtractMatrix(A21, A11), addMatrix(B11, B12));
    auto P7 = strassenMultiply(subtractMatrix(A12, A22), addMatrix(B21, B22));

    auto C11 = addMatrix(subtractMatrix(addMatrix(P1, P4), P5), P7);
    auto C12 = addMatrix(P3, P5);
    auto C21 = addMatrix(P2, P4);
    auto C22 = addMatrix(subtractMatrix(addMatrix(P1, P3), P2), P6);

    vector<vector<T>> C(newSize, vector<T>(newSize));
    for (int i = 0; i < new_dim; ++i)
        for (int j = 0; j < new_dim; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + new_dim] = C12[i][j];
            C[i + new_dim][j] = C21[i][j];
            C[i + new_dim][j + new_dim] = C22[i][j];
        }

    return padded ? unpadMatrix(C, n) : C;
}

template <typename T>
void printMatrix(const vector<vector<T>> &matrix)
{
    int n = matrix.size();
    int m = matrix[0].size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// int
template vector<vector<int>> addMatrix<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> subtractMatrix<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template vector<vector<int>> strassenMultiply<int>(const vector<vector<int>> &A, const vector<vector<int>> &B);
template void printMatrix<int>(const vector<vector<int>> &matrix);
// double
template vector<vector<double>> addMatrix<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> subtractMatrix<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template vector<vector<double>> strassenMultiply<double>(const vector<vector<double>> &A, const vector<vector<double>> &B);
template void printMatrix<double>(const vector<vector<double>> &matrix);
