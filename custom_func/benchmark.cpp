#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>
#include <omp.h>
#include <cblas.h>
#include "custom_trmm.h"

using namespace std;
using namespace std::chrono;

template <typename T>
void fill_random(vector<T>& data) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<T> dist(-1.0, 1.0);
    for (auto& val : data) val = dist(gen);
}

template <typename T>
double run_custom_trmm(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                      CBLAS_DIAG diag, int m, int n, T alpha,
                      const T* A, int lda, const T* B_orig, int ldb, int n_runs) {
    vector<T> B(ldb * n);
    copy(B_orig, B_orig + ldb * n, B.begin());
    
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n_runs; ++i) {
        copy(B_orig, B_orig + ldb * n, B.begin());
        custom_trmm(side, uplo, transA, diag, m, n, alpha, A, lda, B.data(), ldb);
    }
    auto end = high_resolution_clock::now();
    return duration<double>(end - start).count() / n_runs;
}

template <typename T>
double run_openblas_trmm(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                        CBLAS_DIAG diag, int m, int n, T alpha,
                        const T* A, int lda, const T* B_orig, int ldb, int n_runs, int nthreads) {
    openblas_set_num_threads(nthreads);
    vector<T> B(ldb * n);
    copy(B_orig, B_orig + ldb * n, B.begin());
    
    CBLAS_ORDER order = CblasRowMajor;
    auto start = high_resolution_clock::now();
    for (int i = 0; i < n_runs; ++i) {
        copy(B_orig, B_orig + ldb * n, B.begin());
        cblas_dtrmm(order, side, uplo, transA, diag, m, n, alpha, A, lda, B.data(), ldb);
    }
    auto end = high_resolution_clock::now();
    return duration<double>(end - start).count() / n_runs;
}

int main() {
    cout << "TRMM Performance Benchmark" << endl;
    cout << "===========================" << endl;
    
    // Параметры теста
    const int m = 2000, n = 2000;
    const int lda = m, ldb = m;
    const double alpha = 2.0;
    const int n_runs = 10;
    
    cout << "Matrix size: " << m << "x" << n << endl;
    cout << "Number of runs: " << n_runs << endl << endl;
    
    // Создание тестовых данных
    vector<double> A(lda * m), B(ldb * n);
    fill_random(A);
    fill_random(B);
    
    cout << "Testing custom TRMM implementation..." << endl;
    double t_custom = run_custom_trmm(CblasLeft, CblasUpper, CblasNoTrans,
                                     CblasNonUnit, m, n, alpha, A.data(), lda, B.data(), ldb, 1);
    cout << "Custom TRMM time (single run): " << t_custom << " seconds" << endl << endl;
    
    cout << "Thread  | OpenBLAS (s) | Custom (s) | Performance (%)" << endl;
    cout << "--------|--------------|------------|----------------" << endl;
    
    vector<int> threads = {1, 2, 4, 8, 16};
    vector<double> perf_ratios;
    
    for (int nthreads : threads) {
        double t_blas = run_openblas_trmm(CblasLeft, CblasUpper, CblasNoTrans,
                                         CblasNonUnit, m, n, alpha, A.data(), lda, B.data(), ldb, n_runs, nthreads);
        double perf_ratio = (t_custom / t_blas) * 100.0;
        perf_ratios.push_back(perf_ratio);
        
        cout << setw(7) << nthreads << " | " 
             << setw(12) << fixed << setprecision(6) << t_blas << " | "
             << setw(10) << t_custom << " | "
             << setw(15) << setprecision(2) << perf_ratio << endl;
    }
    
    // Вычисление среднего геометрического
    double log_sum = 0;
    for (double ratio : perf_ratios) {
        log_sum += log(ratio);
    }
    double geom_mean = exp(log_sum / perf_ratios.size());
    
    cout << "\nGeometric mean performance: " << geom_mean << "% of OpenBLAS" << endl;
    
    return 0;
}