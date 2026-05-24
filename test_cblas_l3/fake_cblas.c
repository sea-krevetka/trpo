#include <stdio.h>
#include <string.h>
#include "cblas.h"

//Игнорирование параметра транспонирования
void cblas_sgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const float alpha,
                 const float *A, const int lda,
                 const float *B, const int ldb,
                 const float beta,
                 float *C, const int ldc) {
    
    printf("FAKE sgemm called (with bugs!)\n");
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0;
            for (int k = 0; k < K; k++) {
                sum += A[i * lda + k] * B[k * ldb + j];
            }
            C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
        }
    }
}

//Перепутаны M и N
void cblas_dgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const double alpha,
                 const double *A, const int lda,
                 const double *B, const int ldb,
                 const double beta,
                 double *C, const int ldc) {
    
    printf("FAKE dgemm called (M/N swapped!)\n");
    
    //сегфолт или неверный результат
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            double sum = 0;
            for (int k = 0; k < K; k++) {
                if (i < M && k < K && j < N && k < ldb) {
                    sum += A[i * lda + k] * B[k * ldb + j];
                }
            }
            if (i < ldc && j < ldc) {
                C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
            }
        }
    }
}

//Неправильная обработка симметричной матрицы
void cblas_ssymm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const int M, const int N,
                 const float alpha,
                 const float *A, const int lda,
                 const float *B, const int ldb,
                 const float beta,
                 float *C, const int ldc) {
    
    printf("FAKE ssymm called (ignoring symmetry!)\n");
    
    //используем только нижний треугольник
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0;
            if (Side == CblasLeft) {
                // C = A * B
                for (int k = 0; k < M; k++) {
                    sum += A[i * lda + k] * B[k * ldb + j];
                }
            } else {
                // C = B * A
                for (int k = 0; k < N; k++) {
                    sum += B[i * ldb + k] * A[k * lda + j];
                }
            }
            C[i * ldc + j] = alpha * sum + beta * C[i * ldc + j];
        }
    }
}

//Неправильная обработка треугольной матрицы в TRMM
void cblas_strmm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const float alpha,
                 const float *A, const int lda,
                 float *B, const int ldb) {
    
    printf("FAKE strmm called (ignoring triangular structure!)\n");
    
    //используем все элементы
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0;
            if (Side == CblasLeft) {
                // B = A * B
                for (int k = 0; k < M; k++) {
                    sum += A[i * lda + k] * B[k * ldb + j];
                }
            } else {
                // B = B * A
                for (int k = 0; k < N; k++) {
                    sum += B[i * ldb + k] * A[k * lda + j];
                }
            }
            B[i * ldb + j] = alpha * sum;
        }
    }
}

//Неправильное решение в TRSM
void cblas_strsm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const float alpha,
                 const float *A, const int lda,
                 float *B, const int ldb) {
    
    printf("FAKE strsm called (wrong solve!)\n");
    
    
    //копируем A*B вместо решения
    float* temp = (float*)malloc(M * N * sizeof(float));
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0;
            for (int k = 0; k < M; k++) {
                sum += A[i * lda + k] * B[k * ldb + j];
            }
            temp[i * N + j] = sum;
        }
    }
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            B[i * ldb + j] = alpha * temp[i * N + j];
        }
    }
    
    free(temp);
}

//Игнорирование параметра alpha
void cblas_ssyrk(const enum CBLAS_ORDER Order,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans,
                 const int N, const int K,
                 const float alpha,
                 const float *A, const int lda,
                 const float beta,
                 float *C, const int ldc) {
    
    printf("FAKE ssyrk called (ignoring alpha!)\n");
    
    //Всегда alpha=1
    float fake_alpha = 1.0f;
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0;
            for (int k = 0; k < K; k++) {
                if (Trans == CblasNoTrans) {
                    sum += A[i * lda + k] * A[j * lda + k];
                } else {
                    sum += A[k * lda + i] * A[k * lda + j];
                }
            }
            C[i * ldc + j] = fake_alpha * sum + beta * C[i * ldc + j];
        }
    }
}

//и хватит
void cblas_dsymm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const int M, const int N,
                 const double alpha,
                 const double *A, const int lda,
                 const double *B, const int ldb,
                 const double beta,
                 double *C, const int ldc) {
    printf("FAKE dsymm called (not implemented properly)\n");
}

void cblas_cgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const void *alpha,
                 const void *A, const int lda,
                 const void *B, const int ldb,
                 const void *beta,
                 void *C, const int ldc) {
    printf("FAKE cgemm called (complex not implemented)\n");
}

void cblas_zgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const void *alpha,
                 const void *A, const int lda,
                 const void *B, const int ldb,
                 const void *beta,
                 void *C, const int ldc) {
    printf("FAKE zgemm called (complex not implemented)\n");
}

void cblas_dsyrk(const enum CBLAS_ORDER Order,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans,
                 const int N, const int K,
                 const double alpha,
                 const double *A, const int lda,
                 const double beta,
                 double *C, const int ldc) {
    printf("FAKE dsyrk called\n");
}

void cblas_ssyr2k(const enum CBLAS_ORDER Order,
                  const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans,
                  const int N, const int K,
                  const float alpha,
                  const float *A, const int lda,
                  const float *B, const int ldb,
                  const float beta,
                  float *C, const int ldc) {
    printf("FAKE ssyr2k called\n");
}

void cblas_dsyr2k(const enum CBLAS_ORDER Order,
                  const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans,
                  const int N, const int K,
                  const double alpha,
                  const double *A, const int lda,
                  const double *B, const int ldb,
                  const double beta,
                  double *C, const int ldc) {
    printf("FAKE dsyr2k called\n");
}

void cblas_dtrmm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const double alpha,
                 const double *A, const int lda,
                 double *B, const int ldb) {
    printf("FAKE dtrmm called\n");
}

void cblas_dtrsm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const double alpha,
                 const double *A, const int lda,
                 double *B, const int ldb) {
    printf("FAKE dtrsm called\n");
}