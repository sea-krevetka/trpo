#ifndef CBLAS_H
#define CBLAS_H

#include <stddef.h>

enum CBLAS_ORDER {
    CblasRowMajor = 101,
    CblasColMajor = 102
};

enum CBLAS_TRANSPOSE {
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113
};

enum CBLAS_UPLO {
    CblasUpper = 121,
    CblasLower = 122
};

enum CBLAS_DIAG {
    CblasNonUnit = 131,
    CblasUnit = 132
};

enum CBLAS_SIDE {
    CblasLeft = 141,
    CblasRight = 142
};

#ifdef __cplusplus
extern "C" {
#endif

void cblas_sgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const float alpha,
                 const float *A, const int lda,
                 const float *B, const int ldb,
                 const float beta,
                 float *C, const int ldc);

void cblas_dgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const double alpha,
                 const double *A, const int lda,
                 const double *B, const int ldb,
                 const double beta,
                 double *C, const int ldc);

void cblas_cgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const void *alpha,
                 const void *A, const int lda,
                 const void *B, const int ldb,
                 const void *beta,
                 void *C, const int ldc);

void cblas_zgemm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB,
                 const int M, const int N, const int K,
                 const void *alpha,
                 const void *A, const int lda,
                 const void *B, const int ldb,
                 const void *beta,
                 void *C, const int ldc);

// ============= SYMM =============
void cblas_ssymm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const int M, const int N,
                 const float alpha,
                 const float *A, const int lda,
                 const float *B, const int ldb,
                 const float beta,
                 float *C, const int ldc);

void cblas_dsymm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const int M, const int N,
                 const double alpha,
                 const double *A, const int lda,
                 const double *B, const int ldb,
                 const double beta,
                 double *C, const int ldc);

// ============= SYRK =============
void cblas_ssyrk(const enum CBLAS_ORDER Order,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans,
                 const int N, const int K,
                 const float alpha,
                 const float *A, const int lda,
                 const float beta,
                 float *C, const int ldc);

void cblas_dsyrk(const enum CBLAS_ORDER Order,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans,
                 const int N, const int K,
                 const double alpha,
                 const double *A, const int lda,
                 const double beta,
                 double *C, const int ldc);

// ============= SYR2K =============
void cblas_ssyr2k(const enum CBLAS_ORDER Order,
                  const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans,
                  const int N, const int K,
                  const float alpha,
                  const float *A, const int lda,
                  const float *B, const int ldb,
                  const float beta,
                  float *C, const int ldc);

void cblas_dsyr2k(const enum CBLAS_ORDER Order,
                  const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans,
                  const int N, const int K,
                  const double alpha,
                  const double *A, const int lda,
                  const double *B, const int ldb,
                  const double beta,
                  double *C, const int ldc);

// ============= TRMM =============
void cblas_strmm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const float alpha,
                 const float *A, const int lda,
                 float *B, const int ldb);

void cblas_dtrmm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const double alpha,
                 const double *A, const int lda,
                 double *B, const int ldb);

// ============= TRSM =============
void cblas_strsm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const float alpha,
                 const float *A, const int lda,
                 float *B, const int ldb);

void cblas_dtrsm(const enum CBLAS_ORDER Order,
                 const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag,
                 const int M, const int N,
                 const double alpha,
                 const double *A, const int lda,
                 double *B, const int ldb);

#ifdef USE_OPENBLAS
void openblas_set_num_threads(int num_threads);
int openblas_get_num_procs(void);
#endif

#ifdef __cplusplus
}
#endif

#endif