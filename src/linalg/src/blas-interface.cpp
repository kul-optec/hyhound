#include <hyhound/linalg/blas-interface.hpp>
#include <hyhound/linalg/export.h>
#include <hyhound/linalg/lapack.hpp>
#include <hyhound/openmp.h>

namespace hyhound::linalg {

template <>
HYHOUND_LINALG_EXPORT void xgemv(CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA,
                                  blas_index_t M, blas_index_t N, double alpha,
                                  const double *A, blas_index_t lda,
                                  const double *X, blas_index_t incX,
                                  double beta, double *Y, blas_index_t incY) {
    return cblas_dgemv(Layout, TransA, M, N, alpha, A, lda, X, incX, beta, Y,
                       incY);
}
template <>
HYHOUND_LINALG_EXPORT void xgemv(CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA,
                                  blas_index_t M, blas_index_t N, float alpha,
                                  const float *A, blas_index_t lda,
                                  const float *X, blas_index_t incX, float beta,
                                  float *Y, blas_index_t incY) {
    return cblas_sgemv(Layout, TransA, M, N, alpha, A, lda, X, incX, beta, Y,
                       incY);
}

template <>
HYHOUND_LINALG_EXPORT void xgemm(CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA,
                                  CBLAS_TRANSPOSE TransB, index_t M, index_t N,
                                  index_t K, double alpha, const double *A,
                                  index_t lda, const double *B, index_t ldb,
                                  double beta, double *C, index_t ldc) {
    cblas_dgemm(Layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C,
                ldc);
}
template <>
HYHOUND_LINALG_EXPORT void
xgemm(CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
      index_t M, index_t N, index_t K, float alpha, const float *A, index_t lda,
      const float *B, index_t ldb, float beta, float *C, index_t ldc) {
    cblas_sgemm(Layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C,
                ldc);
}

template <>
HYHOUND_LINALG_EXPORT void
xgemmt(CBLAS_LAYOUT Layout, CBLAS_UPLO uplo, CBLAS_TRANSPOSE TransA,
       CBLAS_TRANSPOSE TransB, index_t N, index_t K, double alpha,
       const double *A, index_t lda, const double *B, index_t ldb, double beta,
       double *C, index_t ldc) {
    cblas_dgemmt(Layout, uplo, TransA, TransB, N, K, alpha, A, lda, B, ldb,
                 beta, C, ldc);
}
template <>
HYHOUND_LINALG_EXPORT void
xgemmt(CBLAS_LAYOUT Layout, CBLAS_UPLO uplo, CBLAS_TRANSPOSE TransA,
       CBLAS_TRANSPOSE TransB, index_t N, index_t K, float alpha,
       const float *A, index_t lda, const float *B, index_t ldb, float beta,
       float *C, index_t ldc) {
    cblas_sgemmt(Layout, uplo, TransA, TransB, N, K, alpha, A, lda, B, ldb,
                 beta, C, ldc);
}

template <>
void xtrmv(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
           CBLAS_DIAG Diag, index_t N, const double *A, index_t lda, double *X,
           index_t incX) {
    cblas_dtrmv(Layout, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <>
void xtrmv(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
           CBLAS_DIAG Diag, index_t N, const float *A, index_t lda, float *X,
           index_t incX) {
    cblas_strmv(Layout, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <>
void xtrsv(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
           CBLAS_DIAG Diag, index_t N, const double *A, index_t lda, double *X,
           index_t incX) {
    cblas_dtrsv(Layout, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <>
void xtrsv(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
           CBLAS_DIAG Diag, index_t N, const float *A, index_t lda, float *X,
           index_t incX) {
    cblas_strsv(Layout, Uplo, TransA, Diag, N, A, lda, X, incX);
}

template <>
HYHOUND_LINALG_EXPORT void
xtrmm(CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
      CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N,
      double alpha, const double *A, index_t lda, double *B, index_t ldb) {
    cblas_dtrmm(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}
template <>
HYHOUND_LINALG_EXPORT void
xtrmm(CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
      CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N,
      float alpha, const float *A, index_t lda, float *B, index_t ldb) {
    cblas_strmm(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <>
HYHOUND_LINALG_EXPORT void xsyrk(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo,
                                  CBLAS_TRANSPOSE Trans, index_t N, index_t K,
                                  double alpha, const double *A, index_t lda,
                                  double beta, double *C, index_t ldc) {
    cblas_dsyrk(Layout, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}
template <>
HYHOUND_LINALG_EXPORT void xsyrk(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo,
                                  CBLAS_TRANSPOSE Trans, index_t N, index_t K,
                                  float alpha, const float *A, index_t lda,
                                  float beta, float *C, index_t ldc) {
    cblas_ssyrk(Layout, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}

template <>
HYHOUND_LINALG_EXPORT void
xtrsm(CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
      CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N,
      double alpha, const double *A, index_t lda, double *B, index_t ldb) {
    cblas_dtrsm(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}
template <>
HYHOUND_LINALG_EXPORT void
xtrsm(CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
      CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N,
      float alpha, const float *A, index_t lda, float *B, index_t ldb) {
    cblas_strsm(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}

template <>
void xsytrf_rk(const char *uplo, const index_t *n, double *a,
               const index_t *lda, double *e, index_t *ipiv, double *work,
               const index_t *lwork, index_t *info) {
    dsytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info);
}

template <>
void xsytrf_rk(const char *uplo, const index_t *n, float *a, const index_t *lda,
               float *e, index_t *ipiv, float *work, const index_t *lwork,
               index_t *info) {
    ssytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info);
}

template <>
void xtrtrs(const char *uplo, const char *trans, const char *diag,
            const index_t *n, const index_t *nrhs, const double *A,
            const index_t *ldA, double *B, const index_t *ldB, index_t *info) {
    dtrtrs(uplo, trans, diag, n, nrhs, A, ldA, B, ldB, info);
}

template <>
void xtrtrs(const char *uplo, const char *trans, const char *diag,
            const index_t *n, const index_t *nrhs, const float *A,
            const index_t *ldA, float *B, const index_t *ldB, index_t *info) {
    strtrs(uplo, trans, diag, n, nrhs, A, ldA, B, ldB, info);
}

template <>
void xscal(index_t N, double alpha, double *X, index_t incX) {
    cblas_dscal(N, alpha, X, incX);
}

template <>
void xscal(index_t N, float alpha, float *X, index_t incX) {
    cblas_sscal(N, alpha, X, incX);
}

template <>
void HYHOUND_LINALG_EXPORT xpotrf(const char *uplo, index_t n, double *a,
                                   index_t lda, index_t *info) {
    dpotrf(uplo, &n, a, &lda, info);
}
template <>
void HYHOUND_LINALG_EXPORT xpotrf(const char *uplo, index_t n, float *a,
                                   index_t lda, index_t *info) {
    spotrf(uplo, &n, a, &lda, info);
}

template <>
void HYHOUND_LINALG_EXPORT xlauum(const char *uplo, index_t n, double *a,
                                   index_t lda, index_t *info) {
    dlauum(uplo, &n, a, &lda, info);
}
template <>
void HYHOUND_LINALG_EXPORT xlauum(const char *uplo, index_t n, float *a,
                                   index_t lda, index_t *info) {
    slauum(uplo, &n, a, &lda, info);
}

template <>
void HYHOUND_LINALG_EXPORT xtrtri(const char *uplo, const char *diag,
                                   index_t n, double *a, index_t lda,
                                   index_t *info) {
    dtrtri(uplo, diag, &n, a, &lda, info);
}
template <>
void HYHOUND_LINALG_EXPORT xtrtri(const char *uplo, const char *diag,
                                   index_t n, float *a, index_t lda,
                                   index_t *info) {
    strtri(uplo, diag, &n, a, &lda, info);
}

template <class T, class I>
void xgemv_batch_strided(CBLAS_LAYOUT layout, CBLAS_TRANSPOSE trans, I m, I n,
                         T alpha, const T *a, I lda, I stridea, const T *x,
                         I incx, I stridex, T beta, T *y, I incy, I stridey,
                         I batch_size) {
    HYHOUND_OMP(parallel for)
    for (I i = 0; i < batch_size; ++i) {
        auto offset_a = i * stridea;
        auto offset_x = i * stridex;
        auto offset_y = i * stridey;
        xgemv(layout, trans, m, n, alpha, a + offset_a, lda, x + offset_x, incx,
              beta, y + offset_y, incy);
    }
}

#if HYHOUND_WITH_MKL
template <>
void HYHOUND_LINALG_EXPORT xgemv_batch_strided(
    CBLAS_LAYOUT layout, CBLAS_TRANSPOSE trans, index_t m, index_t n,
    double alpha, const double *a, index_t lda, index_t stridea,
    const double *x, index_t incx, index_t stridex, double beta, double *y,
    index_t incy, index_t stridey, index_t batch_size) {
    if (m == 0 || n == 0)
        return;
    cblas_dgemv_batch_strided(layout, trans, m, n, alpha, a, lda, stridea, x,
                              incx, stridex, beta, y, incy, stridey,
                              batch_size);
}

template <>
void HYHOUND_LINALG_EXPORT xgemv_batch_strided(
    CBLAS_LAYOUT layout, CBLAS_TRANSPOSE trans, index_t m, index_t n,
    float alpha, const float *a, index_t lda, index_t stridea, const float *x,
    index_t incx, index_t stridex, float beta, float *y, index_t incy,
    index_t stridey, index_t batch_size) {
    if (m == 0 || n == 0)
        return;
    cblas_sgemv_batch_strided(layout, trans, m, n, alpha, a, lda, stridea, x,
                              incx, stridex, beta, y, incy, stridey,
                              batch_size);
}
#endif
#if !HYHOUND_WITH_MKL
template void HYHOUND_LINALG_EXPORT xgemv_batch_strided<real_t, index_t>(
    CBLAS_LAYOUT layout, CBLAS_TRANSPOSE trans, index_t m, index_t n,
    real_t alpha, const real_t *a, index_t lda, index_t stridea,
    const real_t *x, index_t incx, index_t stridex, real_t beta, real_t *y,
    index_t incy, index_t stridey, index_t batch_size);
#endif

template <class T, class I>
void xgemm_batch_strided(CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB, I M, I N, I K, T alpha,
                         const T *A, I lda, I stridea, const T *B, I ldb,
                         I strideb, T beta, T *C, I ldc, I stridec,
                         I batch_size) {
    HYHOUND_OMP(parallel for)
    for (I i = 0; i < batch_size; ++i) {
        xgemm(Layout, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C,
              ldc);
        A += stridea;
        B += strideb;
        C += stridec;
    }
}

#if HYHOUND_WITH_MKL
template <>
void HYHOUND_LINALG_EXPORT xgemm_batch_strided(
    CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
    index_t M, index_t N, index_t K, double alpha, const double *A, index_t lda,
    index_t stridea, const double *B, index_t ldb, index_t strideb, double beta,
    double *C, index_t ldc, index_t stridec, index_t batch_size) {
    if (M == 0 || N == 0 || K == 0)
        return;
    cblas_dgemm_batch_strided(Layout, TransA, TransB, M, N, K, alpha, A, lda,
                              stridea, B, ldb, strideb, beta, C, ldc, stridec,
                              batch_size);
}

template <>
void HYHOUND_LINALG_EXPORT xgemm_batch_strided(
    CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
    index_t M, index_t N, index_t K, float alpha, const float *A, index_t lda,
    index_t stridea, const float *B, index_t ldb, index_t strideb, float beta,
    float *C, index_t ldc, index_t stridec, index_t batch_size) {
    if (M == 0 || N == 0 || K == 0)
        return;
    cblas_sgemm_batch_strided(Layout, TransA, TransB, M, N, K, alpha, A, lda,
                              stridea, B, ldb, strideb, beta, C, ldc, stridec,
                              batch_size);
}
#endif
#if !HYHOUND_WITH_MKL
template void HYHOUND_LINALG_EXPORT xgemm_batch_strided<real_t, index_t>(
    CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB,
    index_t M, index_t N, index_t K, real_t alpha, const real_t *A, index_t lda,
    index_t stridea, const real_t *B, index_t ldb, index_t strideb, real_t beta,
    real_t *C, index_t ldc, index_t stridec, index_t batch_size);
#endif

template <class T, class I>
void xsyrk_batch_strided(CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans, I N, I K, T alpha, const T *A,
                         I lda, I stridea, T beta, T *C, I ldc, I stridec,
                         I batch_size) {
    HYHOUND_OMP(parallel for)
    for (I i = 0; i < batch_size; ++i) {
        xsyrk(Layout, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
        A += stridea;
        C += stridec;
    }
}

#if HYHOUND_WITH_MKL
template <>
void HYHOUND_LINALG_EXPORT xsyrk_batch_strided(
    CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, index_t N,
    index_t K, double alpha, const double *A, index_t lda, index_t stridea,
    double beta, double *C, index_t ldc, index_t stridec, index_t batch_size) {
    cblas_dsyrk_batch_strided(Layout, Uplo, Trans, N, K, alpha, A, lda, stridea,
                              beta, C, ldc, stridec, batch_size);
}

template <>
void HYHOUND_LINALG_EXPORT xsyrk_batch_strided(
    CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, index_t N,
    index_t K, float alpha, const float *A, index_t lda, index_t stridea,
    float beta, float *C, index_t ldc, index_t stridec, index_t batch_size) {
    cblas_ssyrk_batch_strided(Layout, Uplo, Trans, N, K, alpha, A, lda, stridea,
                              beta, C, ldc, stridec, batch_size);
}
#endif
#if !HYHOUND_WITH_MKL
template void HYHOUND_LINALG_EXPORT xsyrk_batch_strided<real_t, index_t>(
    CBLAS_LAYOUT Layout, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans, index_t N,
    index_t K, real_t alpha, const real_t *A, index_t lda, index_t stridea,
    real_t beta, real_t *C, index_t ldc, index_t stridec, index_t batch_size);
#endif

template <class T, class I>
void xtrsm_batch_strided(CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, I M, I N,
                         T alpha, const T *A, I lda, I stridea, T *B, I ldb,
                         I strideb, I batch_size) {
    HYHOUND_OMP(parallel for)
    for (I i = 0; i < batch_size; ++i) {
        xtrsm(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
        A += stridea;
        B += strideb;
    }
}

#if HYHOUND_WITH_MKL
template <>
void HYHOUND_LINALG_EXPORT xtrsm_batch_strided(
    CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
    CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N, double alpha,
    const double *A, index_t lda, index_t stridea, double *B, index_t ldb,
    index_t strideb, index_t batch_size) {
    cblas_dtrsm_batch_strided(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A,
                              lda, stridea, B, ldb, strideb, batch_size);
}

template <>
void HYHOUND_LINALG_EXPORT xtrsm_batch_strided(
    CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
    CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N, float alpha,
    const float *A, index_t lda, index_t stridea, float *B, index_t ldb,
    index_t strideb, index_t batch_size) {
    cblas_strsm_batch_strided(Layout, Side, Uplo, TransA, Diag, M, N, alpha, A,
                              lda, stridea, B, ldb, strideb, batch_size);
}
#endif
#if !HYHOUND_WITH_MKL
template void HYHOUND_LINALG_EXPORT xtrsm_batch_strided<real_t, index_t>(
    CBLAS_LAYOUT Layout, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
    CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag, index_t M, index_t N, real_t alpha,
    const real_t *A, index_t lda, index_t stridea, real_t *B, index_t ldb,
    index_t strideb, index_t batch_size);
#endif

template <class T, class I>
void xpotrf_batch_strided(const char *Uplo, I N, T *A, I lda, I stridea,
                          I batch_size) {
    I info_all = 0;
    HYHOUND_OMP(parallel for reduction(+:info_all))
    for (I i = 0; i < batch_size; ++i) {
        I info   = 0;
        I offset = i * stridea;
        xpotrf(Uplo, N, A + offset, lda, &info);
        if (info > 0)
            info = 0; // Ignore factorization failure
        info_all += info;
    }
    // TODO: proper error handling
    lapack_throw_on_err("xpotrf_batch_strided", info_all);
}

template HYHOUND_LINALG_EXPORT void
xpotrf_batch_strided<real_t, index_t>(const char *Uplo, index_t N, real_t *A,
                                      index_t lda, index_t stridea,
                                      index_t batch_size);

} // namespace hyhound::linalg
