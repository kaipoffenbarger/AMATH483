#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <chrono>
#include <limits>
#include <cblas.h>
#include <lapacke.h>
#include <fstream>

using namespace std;

using dcomplex = complex<double>;

double normInf(const dcomplex* A, int n) {
    double max_row_sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            row_sum += abs(A[j * n + i]);  // column-major access
        }
        if (row_sum > max_row_sum) max_row_sum = row_sum;
    }
    return max_row_sum;
}

double norm2(const dcomplex* x, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += norm(x[i]); // norm(x) = |x|^2
    }
    return sqrt(sum);
}

int main() {
    vector<int> sizes;
    for (int n = 16; n <= 8192; n *= 2) sizes.push_back(n);

    ofstream fout("zgesv_errors.txt");
    fout << "# n log10_residual log10_normalized_error\n";

    for (int n : sizes) {
        int ma = n, na = n, nrhs = 1;

        dcomplex* a = (dcomplex*)malloc(sizeof(dcomplex) * ma * na);
        dcomplex* b = (dcomplex*)malloc(sizeof(dcomplex) * ma);
        dcomplex* z = (dcomplex*)malloc(sizeof(dcomplex) * na);
        int* ipiv = (int*)malloc(sizeof(int) * n);

        srand(0);
        int k = 0;
        for (int j = 0; j < na; ++j) {
            for (int i = 0; i < ma; ++i) {
                a[k] = 0.5 - (double)rand() / RAND_MAX
                     + dcomplex(0, 1) * (0.5 - (double)rand() / RAND_MAX);
                if (i == j) a[k] *= static_cast<double>(ma);
                ++k;
            }
        }

        srand(1);
        for (int i = 0; i < ma; ++i) {
            b[i] = 0.5 - (double)rand() / RAND_MAX
                 + dcomplex(0, 1) * (0.5 - (double)rand() / RAND_MAX);
        }

        // Copy original A and b
        vector<dcomplex> A_orig(a, a + n * n);
        vector<dcomplex> b_orig(b, b + n);
        memcpy(z, b, sizeof(dcomplex) * n);

        // Solve A * z = b
        int info = LAPACKE_zgesv(LAPACK_COL_MAJOR, n, nrhs,
                                 reinterpret_cast<lapack_complex_double*>(a), n,
                                 ipiv,
                                 reinterpret_cast<lapack_complex_double*>(z), n);
        if (info != 0) {
            cerr << "LAPACKE_zgesv failed for n = " << n << ", info = " << info << endl;
            continue;
        }

        // Compute Az
        vector<dcomplex> Az(n, dcomplex(0.0, 0.0));
        dcomplex one(1.0, 0.0), zero(0.0, 0.0);
        cblas_zgemv(CblasColMajor, CblasNoTrans, n, n,
                    reinterpret_cast<const void*>(&one), A_orig.data(), n,
                    reinterpret_cast<const void*>(z), 1,
                    reinterpret_cast<const void*>(&zero), Az.data(), 1);

        // Compute residual
        for (int i = 0; i < n; ++i)
            Az[i] -= b_orig[i];
        double residual = norm2(Az.data(), n);

        // Normalized error
        double normA = normInf(A_orig.data(), n);
        double normz = norm2(z, n);
        double eps = numeric_limits<double>::epsilon();
        double normalized_error = residual / (normA * normz * eps);

        fout << n << " "
             << log10(residual) << " "
             << log10(normalized_error) << "\n";

        free(a); free(b); free(z); free(ipiv);
    }

    fout.close();
    cout << "Benchmark complete. Output saved to zgesv_errors.txt\n";
    return 0;
}
