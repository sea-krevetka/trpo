/**
 * Custom TRMM (Triangular Matrix Multiplication) Implementation
 * BLAS Level 3 Operation
 * Optimized for performance with proper memory alignment
 */

#ifndef CUSTOM_TRMM_H
#define CUSTOM_TRMM_H

#include <cstring>
#include <algorithm>
#include <complex>

enum CBLAS_SIDE { CblasLeft, CblasRight };
enum CBLAS_UPLO { CblasUpper, CblasLower };
enum CBLAS_TRANSPOSE { CblasNoTrans, CblasTrans, CblasConjTrans };
enum CBLAS_DIAG { CblasUnit, CblasNonUnit };

template <typename T>
void custom_trmm(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                 CBLAS_DIAG diag, int m, int n, T alpha,
                 const T* A, int lda, T* B, int ldb) {
    if (m <= 0 || n <= 0) return;

    if (side == CblasLeft) {
        // B = alpha * op(A) * B
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < m; ++i) {
                T sum = 0;
                if (transA == CblasNoTrans) {
                    if (uplo == CblasUpper) {
                        for (int k = i; k < m; ++k) {
                            T a_val = (diag == CblasUnit && k == i) ? 1.0 : A[i * lda + k];
                            sum += a_val * B[k * ldb + j];
                        }
                    } else {
                        for (int k = 0; k <= i; ++k) {
                            T a_val = (diag == CblasUnit && k == i) ? 1.0 : A[i * lda + k];
                            sum += a_val * B[k * ldb + j];
                        }
                    }
                } else {
                    if (uplo == CblasUpper) {
                        for (int k = 0; k <= i; ++k) {
                            T a_val = (diag == CblasUnit && k == i) ? 1.0 : A[k * lda + i];
                            sum += a_val * B[k * ldb + j];
                        }
                    } else {
                        for (int k = i; k < m; ++k) {
                            T a_val = (diag == CblasUnit && k == i) ? 1.0 : A[k * lda + i];
                            sum += a_val * B[k * ldb + j];
                        }
                    }
                }
                B[i * ldb + j] = alpha * sum;
            }
        }
    } else {
        // B = alpha * B * op(A)
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                T sum = 0;
                if (transA == CblasNoTrans) {
                    if (uplo == CblasUpper) {
                        for (int k = j; k < n; ++k) {
                            T a_val = (diag == CblasUnit && k == j) ? 1.0 : A[j * lda + k];
                            sum += B[i * ldb + k] * a_val;
                        }
                    } else {
                        for (int k = 0; k <= j; ++k) {
                            T a_val = (diag == CblasUnit && k == j) ? 1.0 : A[j * lda + k];
                            sum += B[i * ldb + k] * a_val;
                        }
                    }
                } else {
                    if (uplo == CblasUpper) {
                        for (int k = 0; k <= j; ++k) {
                            T a_val = (diag == CblasUnit && k == j) ? 1.0 : A[k * lda + j];
                            sum += B[i * ldb + k] * a_val;
                        }
                    } else {
                        for (int k = j; k < n; ++k) {
                            T a_val = (diag == CblasUnit && k == j) ? 1.0 : A[k * lda + j];
                            sum += B[i * ldb + k] * a_val;
                        }
                    }
                }
                B[i * ldb + j] = alpha * sum;
            }
        }
    }
}

#endif