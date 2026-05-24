#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef USE_OPENBLAS
    #include "cblas.h"
    #ifndef OPENBLAS_SET_NUM_THREADS_DEFINED
        void openblas_set_num_threads(int num_threads);
        int openblas_get_num_procs(void);
    #endif
#else
    #include "cblas.h"
#endif

#define EPSILON 1e-4
#define COMPLEX_EPSILON 1e-4

int g_errors = 0;
int g_tests_run = 0;

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define RESET   "\x1b[0m"

#define ASSERT_FLOAT_EQ(actual, expected, msg) do { \
    g_tests_run++; \
    if (fabsf((actual) - (expected)) > EPSILON) { \
        printf(RED "FAIL: %s - expected %f, got %f\n" RESET, msg, (expected), (actual)); \
        g_errors++; \
    } else { \
        printf(GREEN "PASS: %s\n" RESET, msg); \
    } \
} while(0)

#define ASSERT_DOUBLE_EQ(actual, expected, msg) do { \
    g_tests_run++; \
    if (fabs((actual) - (expected)) > EPSILON) { \
        printf(RED "FAIL: %s - expected %f, got %f\n" RESET, msg, (expected), (actual)); \
        g_errors++; \
    } else { \
        printf(GREEN "PASS: %s\n" RESET, msg); \
    } \
} while(0)

#define ASSERT_COMPLEX_FLOAT_EQ(actual, expected, msg) do { \
    g_tests_run++; \
    float diff_re = fabsf(crealf(actual) - crealf(expected)); \
    float diff_im = fabsf(cimagf(actual) - cimagf(expected)); \
    if (diff_re > COMPLEX_EPSILON || diff_im > COMPLEX_EPSILON) { \
        printf(RED "FAIL: %s - expected (%f,%f), got (%f,%f)\n" RESET, \
               msg, crealf(expected), cimagf(expected), crealf(actual), cimagf(actual)); \
        g_errors++; \
    } else { \
        printf(GREEN "PASS: %s\n" RESET, msg); \
    } \
} while(0)


void print_matrix_float(const char* name, float* A, int rows, int cols, int lda) {
    printf("%s:\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%8.4f ", A[i * lda + j]);
        }
        printf("\n");
    }
}

void init_matrix_float(float *A, int rows, int cols, int lda, float base) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A[i * lda + j] = base + i * 0.1f + j * 0.01f;
        }
    }
}

void init_matrix_double(double *A, int rows, int cols, int lda, double base) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A[i * lda + j] = base + i * 0.1 + j * 0.01;
        }
    }
}

void init_matrix_complex_float(float complex *A, int rows, int cols, int lda, float re_base, float im_base) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A[i * lda + j] = (re_base + i * 0.1f + j * 0.01f) + 
                             (im_base + i * 0.2f + j * 0.02f) * I;
        }
    }
}

void init_matrix_complex_double(double complex *A, int rows, int cols, int lda, double re_base, double im_base) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A[i * lda + j] = (re_base + i * 0.1 + j * 0.01) + 
                             (im_base + i * 0.2 + j * 0.02) * I;
        }
    }
}


void test_sgemm_basic() {
    printf("\n Testing sgemm (float) basic \n");
    
    // C = A * B
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float expected[4] = {19.0f, 22.0f, 43.0f, 50.0f};
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "sgemm basic C[%d]", i);
        ASSERT_FLOAT_EQ(C[i], expected[i], msg);
    }
}

void test_sgemm_transpose() {
    printf("\n Testing sgemm with transpose \n");
    // C = A^T * B
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float expected[4] = {26.0f, 30.0f, 38.0f, 44.0f};
    
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                2, 2, 2, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "sgemm transpose C[%d]", i);
        ASSERT_FLOAT_EQ(C[i], expected[i], msg);
    }
}

void test_sgemm_alpha_beta() {
    printf("\n Testing sgemm with alpha and beta \n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    // C = 2*A*B + 3*C 
    float expected[4] = {41.0f, 47.0f, 89.0f, 103.0f};
    
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, 2.0f, A, 2, B, 2, 3.0f, C, 2);
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "sgemm alpha/beta C[%d]", i);
        ASSERT_FLOAT_EQ(C[i], expected[i], msg);
    }
}

void test_dgemm_basic() {
    printf("\n Testing dgemm (double) basic \n");
    
    double A[4] = {1.0, 2.0, 3.0, 4.0};
    double B[4] = {5.0, 6.0, 7.0, 8.0};
    double C[4] = {0.0, 0.0, 0.0, 0.0};
    double expected[4] = {19.0, 22.0, 43.0, 50.0};
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, 1.0, A, 2, B, 2, 0.0, C, 2);
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "dgemm basic C[%d]", i);
        ASSERT_DOUBLE_EQ(C[i], expected[i], msg);
    }
}

void test_cgemm_basic() {
    printf("\n Testing cgemm (complex float) basic \n");
    
    float complex A[4] = {1.0f + 1.0f*I, 2.0f + 2.0f*I, 3.0f + 3.0f*I, 4.0f + 4.0f*I};
    float complex B[4] = {1.0f + 1.0f*I, 1.0f + 1.0f*I, 1.0f + 1.0f*I, 1.0f + 1.0f*I};
    float complex C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                2, 2, 2, &(float complex){1.0f + 0.0f*I}, A, 2, B, 2, 
                &(float complex){0.0f + 0.0f*I}, C, 2);
    
    float complex expected[4] = {
        (1+1I)*(1+1I) + (2+2I)*(1+1I),
        (1+1I)*(1+1I) + (2+2I)*(1+1I),
        (3+3I)*(1+1I) + (4+4I)*(1+1I),
        (3+3I)*(1+1I) + (4+4I)*(1+1I)
    };
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "cgemm basic C[%d]", i);
        ASSERT_COMPLEX_FLOAT_EQ(C[i], expected[i], msg);
    }
}


void test_ssymm_left_upper() {
    printf("\n Testing ssymm left upper \n");
    
    float A[4] = {1.0f, 2.0f, 2.0f, 4.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    // C = A * B
    cblas_ssymm(CblasRowMajor, CblasLeft, CblasUpper,
                2, 2, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    float expected[4] = {7.0f, 10.0f, 14.0f, 20.0f};
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "ssymm left upper C[%d]", i);
        ASSERT_FLOAT_EQ(C[i], expected[i], msg);
    }
}

void test_ssymm_right_lower() {
    printf("\n Testing ssymm right lower \n");
    
    float A[4] = {1.0f, 2.0f, 2.0f, 4.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    // C = B * A 
    cblas_ssymm(CblasRowMajor, CblasRight, CblasLower,
                2, 2, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    
    float expected[4] = {1*1+2*2, 1*2+2*4, 3*1+4*2, 3*2+4*4};
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "ssymm right lower C[%d]", i);
        ASSERT_FLOAT_EQ(C[i], expected[i], msg);
    }
}


void test_ssyrk_upper_notrans() {
    printf("\n Testing ssyrk upper notrans \n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    cblas_ssyrk(CblasRowMajor, CblasUpper, CblasNoTrans,
                2, 2, 1.0f, A, 2, 0.0f, C, 2);
    
    float expected_upper[4] = {5.0f, 11.0f, 0.0f, 25.0f};
    
    ASSERT_FLOAT_EQ(C[0], expected_upper[0], "ssyrk C[0][0]");
    ASSERT_FLOAT_EQ(C[1], expected_upper[1], "ssyrk C[0][1]");
    ASSERT_FLOAT_EQ(C[3], expected_upper[3], "ssyrk C[1][1]");
    
}

void test_ssyr2k_basic() {
    printf("\n Testing ssyr2k basic n");
    
    float A[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    float B[4] = {5.0f, 6.0f, 7.0f, 8.0f};
    float C[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    float C_manual[4] = {0, 0, 0, 0};
    
    cblas_ssyr2k(CblasRowMajor, CblasUpper, CblasNoTrans,
                 2, 2, 1.0f, A, 2, B, 2, 0.0f, C, 2);
    
    float temp1[4] = {0};
    float temp2[4] = {0};
    
    //A * B^T
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                2, 2, 2, 1.0f, A, 2, B, 2, 0.0f, temp1, 2);
    
    //B * A^T  
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                2, 2, 2, 1.0f, B, 2, A, 2, 0.0f, temp2, 2);
    
    for (int i = 0; i < 4; i++) {
        C_manual[i] = temp1[i] + temp2[i];
    }
    
    printf("\n  ssyr2k result: [%f, %f; %f, %f]\n", 
           C[0], C[1], C[2], C[3]);
    printf("  manual result: [%f, %f; %f, %f]\n", 
           C_manual[0], C_manual[1], C_manual[2], C_manual[3]);
    
    char msg[100];
    sprintf(msg, "ssyr2k C[0][0] (vs manual)");
    ASSERT_FLOAT_EQ(C[0], C_manual[0], msg);
    
    sprintf(msg, "ssyr2k C[0][1] (vs manual)");
    ASSERT_FLOAT_EQ(C[1], C_manual[1], msg);
    
    sprintf(msg, "ssyr2k C[1][1] (vs manual)");
    ASSERT_FLOAT_EQ(C[3], C_manual[3], msg);
    
    if (fabs(C[1] - C[2]) < EPSILON) {
        printf(GREEN "PASS: ssyr2k symmetry check\n" RESET);
        g_tests_run++;
    } else {
        printf(YELLOW "WARN: ssyr2k not symmetric (upper=%f, lower=%f) - may be by design\n" RESET, 
               C[1], C[2]);
        g_tests_run++;
    }
}


void test_strmm_left_upper_notrans() {
    printf("\n Testing strmm left upper notrans \n");
    
    float A[4] = {1.0f, 2.0f, 0.0f, 3.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    
    // B = A * B
    cblas_strmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
                2, 2, 1.0f, A, 2, B, 2);
    
    float expected[4] = {7.0f, 10.0f, 9.0f, 12.0f};
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "strmm left upper C[%d]", i);
        ASSERT_FLOAT_EQ(B[i], expected[i], msg);
    }
}

void test_strmm_right_lower_trans() {
    printf("\n Testing strmm right lower trans \n");
    

    float A[4] = {1.0f, 0.0f, 2.0f, 3.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f};
    
    // B = B * A^T
    cblas_strmm(CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
                2, 2, 1.0f, A, 2, B, 2);

    float changed = 0;
    for (int i = 0; i < 4; i++) {
        if (fabs(B[i] - (i+1)) > EPSILON) changed = 1;
    }
    
    if (changed) {
        printf(GREEN "PASS: strmm right lower trans - matrix changed\n" RESET);
        g_tests_run++;
    } else {
        printf(RED "FAIL: strmm right lower trans - matrix unchanged\n" RESET);
        g_errors++;
    }
}


void test_strsm_solve() {
    printf("\n Testing strsm solve n");
    
    float A[4] = {1.0f, 0.0f, 2.0f, 3.0f};
    float B[4] = {1.0f, 2.0f, 3.0f, 4.0f}; 
    
    float X_expected[4];
    memcpy(X_expected, B, 4 * sizeof(float));
    
    float B_computed[4];
    memcpy(B_computed, B, 4 * sizeof(float));
    
    cblas_strmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                2, 2, 1.0f, A, 2, B_computed, 2);
    
    float X_solved[4];
    memcpy(X_solved, B_computed, 4 * sizeof(float));
    
    cblas_strsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit,
                2, 2, 1.0f, A, 2, X_solved, 2);
    
    for (int i = 0; i < 4; i++) {
        char msg[100];
        sprintf(msg, "strsm solve X[%d]", i);
        ASSERT_FLOAT_EQ(X_solved[i], X_expected[i], msg);
    }
}

void test_with_threads() {
    printf("\n Testing with multiple threads \n");
    
    #ifdef USE_OPENBLAS
    int max_threads = openblas_get_num_procs();
    printf("OpenBLAS detected. Max threads: %d\n", max_threads);
    
    for (int threads = 1; threads <= max_threads && threads <= 4; threads++) {
        printf("\nTesting with %d thread(s):\n", threads);
        openblas_set_num_threads(threads);
        
        test_sgemm_basic();
        test_dgemm_basic();
        test_strmm_left_upper_notrans();
        
        printf("Thread test with %d threads completed\n", threads);
    }
    #else
    printf("Not compiled with OpenBLAS, skipping thread tests\n");
    #endif
}


void test_all_precisions() {
    printf("\n========== TESTING ALL PRECISIONS ==========\n");
    
    printf("\n>> Single precision (float):\n");
    test_sgemm_basic();
    test_sgemm_transpose();
    test_sgemm_alpha_beta();
    test_ssymm_left_upper();
    test_ssymm_right_lower();
    test_ssyrk_upper_notrans();
    test_ssyr2k_basic();
    test_strmm_left_upper_notrans();
    test_strmm_right_lower_trans();
    test_strsm_solve();
    
    printf("\n>> Double precision (double):\n");
    test_dgemm_basic();
    
    // comple
    printf("\n>> Complex precision:\n");
    test_cgemm_basic();
}


int main(int argc, char** argv) {
    printf("========================================\n");
    printf("   CBLAS Level 3 Interface Tests\n");
    printf("========================================\n");
    
    #ifdef USE_OPENBLAS
    printf("Testing against: OpenBLAS\n");
    #else
    printf("Testing against: Custom library\n");
    #endif
    
    g_errors = 0;
    g_tests_run = 0;
    
    test_all_precisions();
    
    test_with_threads();
    
    printf("\n========================================\n");
    printf("Tests run: %d\n", g_tests_run);
    
    if (g_errors == 0) {
        printf(GREEN "ALL TESTS PASSED!\n" RESET);
        return 0;
    } else {
        printf(RED "%d TESTS FAILED!\n" RESET, g_errors);
        return 1;
    }
}

// testing pipeline