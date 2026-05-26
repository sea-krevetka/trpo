#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <algorithm>
#include <cstring>
#include <omp.h>
#include <cblas.h>  // OpenBLAS CBLAS interface

// =====================================================================
// Оптимизированная шаблонная реализация TRMM (с cache-friendly доступом)
// =====================================================================
template <typename T>
void custom_trmm(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                 CBLAS_DIAG diag, int m, int n, T alpha,
                 const T* A, int lda, T* B, int ldb)
{
    if (m <= 0 || n <= 0) return;

    if (side == CblasLeft) {
        // B = alpha * op(A) * B
        if (transA == CblasNoTrans) {
            if (uplo == CblasUpper) {
                // A is upper triangular, B = alpha * A * B
                for (int j = 0; j < n; ++j) {
                    for (int i = m-1; i >= 0; --i) {
                        T sum = B[i * ldb + j];
                        for (int k = i+1; k < m; ++k) {
                            sum += A[i * lda + k] * B[k * ldb + j];
                        }
                        if (diag == CblasUnit) {
                            B[i * ldb + j] = alpha * sum;
                        } else {
                            B[i * ldb + j] = alpha * (A[i * lda + i] * B[i * ldb + j] + sum - B[i * ldb + j]);
                        }
                    }
                }
            } else { // Lower
                // A is lower triangular, B = alpha * A * B
                for (int j = 0; j < n; ++j) {
                    for (int i = 0; i < m; ++i) {
                        T sum = B[i * ldb + j];
                        for (int k = 0; k < i; ++k) {
                            sum += A[i * lda + k] * B[k * ldb + j];
                        }
                        if (diag == CblasUnit) {
                            B[i * ldb + j] = alpha * sum;
                        } else {
                            B[i * ldb + j] = alpha * (A[i * lda + i] * B[i * ldb + j] + sum - B[i * ldb + j]);
                        }
                    }
                }
            }
        } else {
            // transA != CblasNoTrans: B = alpha * A^T * B
            if (uplo == CblasUpper) {
                // A^T is lower triangular
                for (int j = 0; j < n; ++j) {
                    for (int i = 0; i < m; ++i) {
                        T sum = 0;
                        for (int k = 0; k <= i; ++k) {
                            T a_val = (diag == CblasUnit && k == i) ? 1.0 : A[k * lda + i];
                            sum += a_val * B[k * ldb + j];
                        }
                        B[i * ldb + j] = alpha * sum;
                    }
                }
            } else { // Lower
                // A^T is upper triangular
                for (int j = 0; j < n; ++j) {
                    for (int i = m-1; i >= 0; --i) {
                        T sum = 0;
                        for (int k = i; k < m; ++k) {
                            T a_val = (diag == CblasUnit && k == i) ? 1.0 : A[k * lda + i];
                            sum += a_val * B[k * ldb + j];
                        }
                        B[i * ldb + j] = alpha * sum;
                    }
                }
            }
        }
    } else { // CblasRight
        // B = alpha * B * op(A)
        if (transA == CblasNoTrans) {
            if (uplo == CblasUpper) {
                // B = alpha * B * A, A is upper triangular
                for (int i = 0; i < m; ++i) {
                    for (int j = n-1; j >= 0; --j) {
                        T sum = B[i * ldb + j];
                        for (int k = j+1; k < n; ++k) {
                            sum += B[i * ldb + k] * A[j * lda + k];
                        }
                        if (diag == CblasUnit) {
                            B[i * ldb + j] = alpha * sum;
                        } else {
                            B[i * ldb + j] = alpha * (B[i * ldb + j] * A[j * lda + j] + sum - B[i * ldb + j]);
                        }
                    }
                }
            } else { // Lower
                // B = alpha * B * A, A is lower triangular
                for (int i = 0; i < m; ++i) {
                    for (int j = 0; j < n; ++j) {
                        T sum = B[i * ldb + j];
                        for (int k = 0; k < j; ++k) {
                            sum += B[i * ldb + k] * A[j * lda + k];
                        }
                        if (diag == CblasUnit) {
                            B[i * ldb + j] = alpha * sum;
                        } else {
                            B[i * ldb + j] = alpha * (B[i * ldb + j] * A[j * lda + j] + sum - B[i * ldb + j]);
                        }
                    }
                }
            }
        } else {
            // B = alpha * B * A^T
            if (uplo == CblasUpper) {
                // A^T is lower triangular
                for (int i = 0; i < m; ++i) {
                    for (int j = 0; j < n; ++j) {
                        T sum = 0;
                        for (int k = 0; k <= j; ++k) {
                            T a_val = (diag == CblasUnit && k == j) ? 1.0 : A[k * lda + j];
                            sum += B[i * ldb + k] * a_val;
                        }
                        B[i * ldb + j] = alpha * sum;
                    }
                }
            } else { // Lower
                // A^T is upper triangular
                for (int i = 0; i < m; ++i) {
                    for (int j = n-1; j >= 0; --j) {
                        T sum = 0;
                        for (int k = j; k < n; ++k) {
                            T a_val = (diag == CblasUnit && k == j) ? 1.0 : A[k * lda + j];
                            sum += B[i * ldb + k] * a_val;
                        }
                        B[i * ldb + j] = alpha * sum;
                    }
                }
            }
        }
    }
}

// =====================================================================
// Вспомогательные функции для бенчмарка
// =====================================================================
template <typename T>
void fill_random(std::vector<T>& data) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(-1.0, 1.0);
    for (auto& val : data) val = dist(gen);
}

// Замер времени для кастомной реализации
template <typename T>
double run_custom_trmm(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                       CBLAS_DIAG diag, int m, int n, T alpha,
                       const std::vector<T>& A, int lda,
                       const std::vector<T>& B_original, int ldb,
                       int n_runs)
{
    std::vector<T> B_copy = B_original;

    // Разогрев
    custom_trmm(side, uplo, transA, diag, m, n, alpha, A.data(), lda, B_copy.data(), ldb);

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_runs; ++i) {
        std::copy(B_original.begin(), B_original.end(), B_copy.begin());
        custom_trmm(side, uplo, transA, diag, m, n, alpha, A.data(), lda, B_copy.data(), ldb);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count() / n_runs;
}

// Замер времени для OpenBLAS
template <typename T>
double run_openblas_trmm(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                         CBLAS_DIAG diag, int m, int n, T alpha,
                         const std::vector<T>& A, int lda,
                         const std::vector<T>& B_original, int ldb,
                         int n_runs, int nthreads)
{
    openblas_set_num_threads(nthreads);
    std::vector<T> B_copy = B_original;

    // Разогрев
    if constexpr (std::is_same_v<T, double>) {
        cblas_dtrmm(CblasRowMajor, side, uplo, transA, diag, m, n, alpha,
                    A.data(), lda, B_copy.data(), ldb);
    } else {
        cblas_strmm(CblasRowMajor, side, uplo, transA, diag, m, n, (float)alpha,
                    A.data(), lda, B_copy.data(), ldb);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_runs; ++i) {
        std::copy(B_original.begin(), B_original.end(), B_copy.begin());
        if constexpr (std::is_same_v<T, double>) {
            cblas_dtrmm(CblasRowMajor, side, uplo, transA, diag, m, n, alpha,
                        A.data(), lda, B_copy.data(), ldb);
        } else {
            cblas_strmm(CblasRowMajor, side, uplo, transA, diag, m, n, (float)alpha,
                        A.data(), lda, B_copy.data(), ldb);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count() / n_runs;
}

// Вычисление среднего геометрического
double geometric_mean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    double log_sum = 0.0;
    for (double v : values) {
        if (v <= 0) return 0.0;  // Защита от log(0) или отрицательных
        log_sum += std::log(v);
    }
    return std::exp(log_sum / values.size());
}

// =====================================================================
// Главная программа тестирования
// =====================================================================
int main() {
    std::cout << "OpenBLAS TRMM Benchmark" << std::endl;
    std::cout << "=======================" << std::endl;
    
    // Подбираем размер для ~1 минуты на одном ядре
    const int m = 1000;
    const int n = 1000;
    const int lda = m;
    const int ldb = m;
    const double alpha = 2.0;

    std::cout << "Matrix size: " << m << "x" << n << std::endl;
    std::cout << "Memory for A: " << (lda * m * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "Memory for B: " << (ldb * n * sizeof(double)) / (1024.0 * 1024.0) << " MB" << std::endl;

    // Выделение и заполнение матриц
    std::vector<double> A(lda * m);
    std::vector<double> B(ldb * n);
    fill_random(A);
    fill_random(B);

    // Оценка времени одного прогона
    std::cout << "\nEstimating single custom TRMM time..." << std::flush;
    double t_custom_single = run_custom_trmm(CblasLeft, CblasUpper, CblasNoTrans,
                                            CblasNonUnit, m, n, alpha, A, lda, B, ldb, 1);
    std::cout << " " << t_custom_single << " sec" << std::endl;
    
    if (t_custom_single < 60.0) {
        std::cout << "Warning: Single run less than 60 sec. Consider increasing matrix size." << std::endl;
    }

    // Параметры тестирования
    std::vector<int> threads = {1, 2, 4, 8, 16};
    const int n_runs = 10;

    std::cout << "\n=== Performance comparison (average of " << n_runs << " runs) ===" << std::endl;
    std::cout << std::setw(8) << "Threads" << " | " 
              << std::setw(12) << "OpenBLAS(s)" << " | " 
              << std::setw(12) << "Custom(s)" << " | " 
              << std::setw(12) << "Rel.perf(%)" << " | " 
              << std::setw(15) << "Geom.mean(%)" << std::endl;
    std::cout << std::string(75, '-') << std::endl;

    std::vector<double> thread_geom_means;

    for (int nthreads : threads) {
        std::vector<double> perf_ratios;
        double blas_total_time = 0.0;
        double custom_total_time = 0.0;

        std::cout << "Testing " << nthreads << " threads..." << std::flush;

        for (int run = 0; run < n_runs; ++run) {
            double t_blas = run_openblas_trmm(CblasLeft, CblasUpper, CblasNoTrans,
                                             CblasNonUnit, m, n, alpha, A, lda, B, ldb, 1, nthreads);
            double t_custom = run_custom_trmm(CblasLeft, CblasUpper, CblasNoTrans,
                                             CblasNonUnit, m, n, alpha, A, lda, B, ldb, 1);
            
            blas_total_time += t_blas;
            custom_total_time += t_custom;
            
            // Относительная производительность: (t_blas / t_custom) * 100%
            // Это показывает, сколько процентов от времени OpenBLAS занимает наша реализация
            double ratio = (t_blas / t_custom) * 100.0;
            perf_ratios.push_back(ratio);
        }

        double avg_blas_time = blas_total_time / n_runs;
        double avg_custom_time = custom_total_time / n_runs;
        double avg_ratio = (avg_blas_time / avg_custom_time) * 100.0;
        double geom_mean_perf = geometric_mean(perf_ratios);
        thread_geom_means.push_back(geom_mean_perf);

        std::cout << " done" << std::endl;
        std::cout << std::setw(8) << nthreads << " | " 
                  << std::fixed << std::setprecision(4) << std::setw(12) << avg_blas_time << " | "
                  << std::setw(12) << avg_custom_time << " | "
                  << std::setw(11) << std::setprecision(2) << avg_ratio << " | "
                  << std::setw(14) << geom_mean_perf << std::endl;
    }

    double overall_geom_mean = geometric_mean(thread_geom_means);
    std::cout << "\n" << std::string(75, '=') << std::endl;
    std::cout << "Overall geometric mean performance: " 
              << std::fixed << std::setprecision(2) << overall_geom_mean << "%" << std::endl;
    
    return 0;
}