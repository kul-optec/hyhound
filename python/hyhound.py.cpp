#include <hyhound/householder-updowndate.hpp>
#include <hyhound/updown.hpp>
#include <hyhound-version.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/tuple.h>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace nb = nanobind;
using namespace nb::literals;

namespace hyhound::py {
template <class T>
constexpr hyhound::Config<T> config{};
template <class T, class UpDown>
void update_cholesky(MatrixView<T> L, MatrixView<T> A, const UpDown &signs,
                     MatrixView<T> Ws = MatrixView<T>{{.rows = 0}}) {
    nb::gil_scoped_release release;
    hyhound::update_cholesky<T, config<T>>(L, A, signs, Ws);
}
template <class T, class UpDown>
void apply_householder(MatrixView<T> L, MatrixView<T> A, const UpDown &signs,
                       std::type_identity_t<MatrixView<const T>> Ws,
                       std::type_identity_t<MatrixView<const T>> B) {
    nb::gil_scoped_release release;
    hyhound::apply_householder<T, config<T>>(L, A, signs, Ws, B);
}
} // namespace hyhound::py

namespace {

template <class... Args>
auto view(const nb::ndarray<Args...> &array)
    requires(array.Order == 'F')
{
    return hyhound::MatrixView<typename nb::ndarray<Args...>::Scalar>{{
        .data         = array.data(),
        .rows         = static_cast<hyhound::index_t>(array.shape(0)),
        .cols         = static_cast<hyhound::index_t>(array.shape(1)),
        .outer_stride = static_cast<hyhound::index_t>(array.stride(1)),
    }};
}

template <class... Args>
auto copy(const nb::ndarray<Args...> &array) {
    using T     = std::remove_cv_t<typename nb::ndarray<Args...>::Scalar>;
    auto v      = array.view();
    size_t rows = v.shape(0), cols = v.shape(1);
    std::unique_ptr<T[]> data{new T[rows * cols]};
    auto deleter = [](void *p) noexcept { delete[] static_cast<T *>(p); };
    if (v.stride(0) <= v.stride(1))
        for (size_t c = 0; c < cols; ++c)
            for (size_t r = 0; r < rows; ++r)
                data[r + c * rows] = v(r, c);
    else
        for (size_t r = 0; r < rows; ++r)
            for (size_t c = 0; c < cols; ++c)
                data[r + c * rows] = v(r, c);
    auto data_ptr = data.get();
    return nb::ndarray<nb::numpy, T, nb::ndim<2>, nb::f_contig>{
        data_ptr,
        {rows, cols},
        nb::capsule{data.release(), deleter},
        {1, static_cast<int64_t>(rows)},
    };
}

template <class T>
auto empty(size_t rows, size_t cols) {
    std::unique_ptr<T[]> data{new T[rows * cols]()};
    auto deleter  = [](void *p) noexcept { delete[] static_cast<T *>(p); };
    auto data_ptr = data.get();
    return nb::ndarray<nb::numpy, T, nb::ndim<2>, nb::f_contig>{
        data_ptr,
        {rows, cols},
        nb::capsule{data.release(), deleter},
        {1, static_cast<int64_t>(rows)},
    };
}

template <class... Args>
auto zero_top_rows(nb::ndarray<Args...> &array, size_t n) {
    using T     = typename nb::ndarray<Args...>::Scalar;
    auto v      = array.view();
    size_t rows = v.shape(0), cols = v.shape(1);
    assert(n <= rows);
    for (size_t c = 0; c < cols; ++c)
        for (size_t r = 0; r < n; ++r)
            v(r, c) = T{};
}

template <class... ArgsL, class... ArgsA>
void check_dim(const nb::ndarray<ArgsL...> &L, const nb::ndarray<ArgsA...> &A) {
    if (L.ndim() != 2)
        throw std::invalid_argument("L.ndim should be 2");
    if (A.ndim() != 2)
        throw std::invalid_argument("A.ndim should be 2");
    if (L.shape(0) < L.shape(1))
        throw std::invalid_argument("L.shape[0] should be greater than "
                                    "or equal to L.shape[1]");
    if (A.shape(0) != L.shape(0))
        throw std::invalid_argument("A.shape[0] should match L.shape[0]");
}

template <class... ArgsL, class... ArgsA, class... ArgsW, class... ArgsB>
void check_dim(const nb::ndarray<ArgsL...> &L, const nb::ndarray<ArgsA...> &A,
               const nb::ndarray<ArgsW...> &W, const nb::ndarray<ArgsB...> &B) {
    check_dim(L, A);
    using T = std::remove_const_t<typename nb::ndarray<ArgsW...>::Scalar>;
    static constexpr auto R = hyhound::py::config<T>.block_size_r;
    if (W.shape(0) != R)
        throw std::invalid_argument("W.shape[0] should be " +
                                    std::to_string(R));
    if (W.shape(1) != L.shape(1))
        throw std::invalid_argument("W.shape[1] should match L.shape[1]");
    if (B.shape(0) != L.shape(1))
        throw std::invalid_argument("B.shape[0] should match L.shape[1]");
    if (B.shape(1) != A.shape(1))
        throw std::invalid_argument("B.shape[1] should match A.shape[1]");
}

template <class T>
void register_module(nb::module_ &m) {
    using c_matrix = nb::ndarray<const T, nb::ndim<2>, nb::device::cpu>;
    using c_vector =
        nb::ndarray<const T, nb::ndim<1>, nb::device::cpu, nb::any_contig>;
    using matrix_cm =
        nb::ndarray<T, nb::ndim<2>, nb::device::cpu, nb::f_contig>;
    using c_matrix_cm =
        nb::ndarray<const T, nb::ndim<2>, nb::device::cpu, nb::f_contig>;
    // In-place
    m.def(
        "update_cholesky_inplace",
        [](matrix_cm L, matrix_cm A) {
            check_dim(L, A);
            hyhound::py::update_cholesky(view(L), view(A), hyhound::Update{});
            zero_top_rows(A, L.shape(1));
        },
        "L"_a.noconvert(), "A"_a.noconvert(),
        R"doc(
Cholesky factorization update. Overwrites its arguments.

L̃L̃ᵀ + ÃÃᵀ = LLᵀ + AAᵀ

Parameters
----------
L : (k × n), lower-trapezoidal, Fortran order
    On entry, the original Cholesky factor L.
    On exit, contains the updated Cholesky factor L̃.

A : (k × m), rectangular, Fortran order
    On entry, the update matrix A.
    On exit, contains the k-n bottom rows of the remaining update matrix Ã
    (the top n rows of Ã are zero).
)doc");

    m.def(
        "downdate_cholesky_inplace",
        [](matrix_cm L, matrix_cm A) {
            check_dim(L, A);
            hyhound::py::update_cholesky(view(L), view(A), hyhound::Downdate{});
            zero_top_rows(A, L.shape(1));
        },
        "L"_a.noconvert(), "A"_a.noconvert(),
        R"doc(
Cholesky factorization downdate. Overwrites its arguments.

L̃L̃ᵀ - ÃÃᵀ = LLᵀ - AAᵀ

Parameters
----------
L : (k × n), lower-trapezoidal, Fortran order
    On entry, the original Cholesky factor L.
    On exit, contains the updated Cholesky factor L̃.

A : (k × m), rectangular, Fortran order
    On entry, the downdate matrix A.
    On exit, contains the k-n bottom rows of the remaining downdate matrix Ã
    (the top n rows of Ã are zero).
)doc");

    m.def(
        "update_cholesky_sign_inplace",
        [](matrix_cm L, matrix_cm A, c_vector signs) {
            check_dim(L, A);
            if (A.shape(1) != signs.size())
                throw std::invalid_argument("len(signs) should be A.shape[1]");
            std::span signs_span{signs.data(), signs.shape(0)};
            if (!std::ranges::all_of(signs_span, [](T t) { return t == T{}; }))
                throw std::invalid_argument("signs should be +/- zero");
            hyhound::UpDowndate<T> sgn{signs_span};
            hyhound::py::update_cholesky(view(L), view(A), sgn);
            zero_top_rows(A, L.shape(1));
        },
        "L"_a.noconvert(), "A"_a.noconvert(), "signs"_a,
        R"doc(
Cholesky factorization update with signed columns. Overwrites its arguments.

L̃L̃ᵀ + ÃSÃᵀ = LLᵀ + ASAᵀ,
where S = np.diag(np.copysign(np.ones(m), signs)) and signs contains ±0.

Parameters
----------
L : (k × n), lower-trapezoidal, Fortran order
    On entry, the original Cholesky factor L.
    On exit, contains the updated Cholesky factor L̃.

A : (k × m), rectangular, Fortran order
    On entry, the update matrix A.
    On exit, contains the k-n bottom rows of the remaining update matrix Ã
    (the top n rows of Ã are zero).

signs : m-vector
    Signs that determine whether a column of A is added (+0) or removed (-0).
    Values other than ±0 are not allowed.
)doc");

    m.def(
        "update_cholesky_diag_inplace",
        [](matrix_cm L, matrix_cm A, c_vector diag) {
            check_dim(L, A);
            if (A.shape(1) != diag.size())
                throw std::invalid_argument("len(diag) should be A.shape[1]");
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::py::update_cholesky(view(L), view(A), d);
            zero_top_rows(A, L.shape(1));
        },
        "L"_a.noconvert(), "A"_a.noconvert(), "diag"_a,
        R"doc(
Cholesky factorization update with diagonal scaling. Overwrites its arguments.

L̃L̃ᵀ + ÃDÃᵀ = LLᵀ + ADAᵀ,
where D = np.diag(diag).

Parameters
----------
L : (k × n), lower-trapezoidal, Fortran order
    On entry, the original Cholesky factor L.
    On exit, contains the updated Cholesky factor L̃.

A : (k × m), rectangular, Fortran order
    On entry, the update matrix A.
    On exit, contains the k-n bottom rows of the remaining update matrix Ã
    (the top n rows of Ã are zero).

diag : m-vector
    Scale factors corresponding to the columns of A.
)doc");

    // Returning copies
    m.def(
        "update_cholesky",
        [](c_matrix L, c_matrix A) {
            check_dim(L, A);
            auto L̃ = copy(L), Ã = copy(A);
            auto W = empty<T>(hyhound::py::config<T>.block_size_r, L̃.shape(1));
            hyhound::py::update_cholesky(view(L̃), view(Ã), hyhound::Update{},
                                         view(W));
            return std::make_tuple(std::move(L̃), std::move(Ã), std::move(W));
        },
        "L"_a, "A"_a,
        R"doc(
Cholesky factorization update. Returns updated copies.

L̃L̃ᵀ + ÃÃᵀ = LLᵀ + AAᵀ

Parameters
----------
L : (k × n), lower-trapezoidal
    The original Cholesky factor.

A : (k × m), rectangular
    The update matrix.

Returns
-------
L̃ : (k × n)
    The updated Cholesky factor.

A_rem : (k × m)
    Contains the k-n bottom rows of the remaining update matrix Ã.
    The top n rows of Ã are zero (not stored explicitly).
    The top n rows of A_rem contain Householder reflectors.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update. Together with the top n rows of A_rem,
    this can be used to apply the block Householder transformation to other matrices.
    The number of rows depends on the block size and is architecture-dependent.
)doc");

    m.def(
        "downdate_cholesky",
        [](c_matrix L, c_matrix A) {
            check_dim(L, A);
            auto L̃ = copy(L), Ã = copy(A);
            auto W = empty<T>(hyhound::py::config<T>.block_size_r, L̃.shape(1));
            hyhound::py::update_cholesky(view(L̃), view(Ã), hyhound::Downdate{},
                                         view(W));
            return std::make_tuple(std::move(L̃), std::move(Ã), std::move(W));
        },
        "L"_a, "A"_a,
        R"doc(
Cholesky factorization downdate. Returns updated copies.

L̃L̃ᵀ - ÃÃᵀ = LLᵀ - AAᵀ

Parameters
----------
L : (k × n), lower-trapezoidal
    The original Cholesky factor.

A : (k × m), rectangular
    The downdate matrix.

Returns
-------
L̃ : (k × n)
    The updated Cholesky factor.

A_rem : (k × m)
    Contains the k-n bottom rows of the remaining downdate matrix Ã.
    The top n rows of Ã are zero (not stored explicitly).
    The top n rows of A_rem contain Householder reflectors.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update. Together with the top n rows of A_rem,
    this can be used to apply the block Householder transformation to other matrices.
    The number of rows depends on the block size and is architecture-dependent.
)doc");

    m.def(
        "update_cholesky_sign",
        [](c_matrix L, c_matrix A, c_vector signs) {
            check_dim(L, A);
            if (A.shape(1) != signs.size())
                throw std::invalid_argument("len(signs) should be A.shape[1]");
            std::span signs_span{signs.data(), signs.shape(0)};
            if (!std::ranges::all_of(signs_span, [](T t) { return t == T{}; }))
                throw std::invalid_argument("signs should be +/- zero");
            auto L̃ = copy(L), Ã = copy(A);
            auto W = empty<T>(hyhound::py::config<T>.block_size_r, L̃.shape(1));
            hyhound::UpDowndate<T> sgn{signs_span};
            hyhound::py::update_cholesky(view(L̃), view(Ã), sgn, view(W));
            return std::make_tuple(std::move(L̃), std::move(Ã), std::move(W));
        },
        "L"_a, "A"_a, "signs"_a,
        R"doc(
Cholesky factorization update with signed columns. Returns updated copies.

L̃L̃ᵀ + ÃSÃᵀ = LLᵀ + ASAᵀ,
where S = np.diag(np.copysign(np.ones(m), signs)) and signs contains ±0.

Parameters
----------
L : (k × n), lower-trapezoidal
    The original Cholesky factor.

A : (k × m), rectangular
    The update matrix.

signs : m-vector
    Signs that determine whether a column of A is added (+0) or removed (-0).
    Values other than ±0 are not allowed.

Returns
-------
L̃ : (k × n)
    The updated Cholesky factor.

A_rem : (k × m)
    Contains the k-n bottom rows of the remaining update matrix Ã.
    The top n rows of Ã are zero (not stored explicitly).
    The top n rows of A_rem contain Householder reflectors.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update. Together with the top n rows of A_rem,
    this can be used to apply the block Householder transformation to other matrices.
    The number of rows depends on the block size and is architecture-dependent.
)doc");

    m.def(
        "update_cholesky_diag",
        [](c_matrix L, c_matrix A, c_vector diag) {
            check_dim(L, A);
            if (A.shape(1) != diag.size())
                throw std::invalid_argument("len(diag) should be A.shape[1]");
            auto L̃ = copy(L), Ã = copy(A);
            auto W = empty<T>(hyhound::py::config<T>.block_size_r, L̃.shape(1));
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::py::update_cholesky(view(L̃), view(Ã), d, view(W));
            return std::make_tuple(std::move(L̃), std::move(Ã), std::move(W));
        },
        "L"_a, "A"_a, "diag"_a,
        R"doc(
Cholesky factorization update with diagonal scaling. Returns updated copies.

L̃L̃ᵀ + ÃDÃᵀ = LLᵀ + ADAᵀ,
where D = np.diag(diag).

Parameters
----------
L : (k × n), lower-trapezoidal
    The original Cholesky factor.

A : (k × m), rectangular
    The update matrix.

diag : m-vector
    Scale factors corresponding to the columns of A.

Returns
-------
L̃ : (k × n)
    The updated Cholesky factor.

A_rem : (k × m)
    Contains the k-n bottom rows of the remaining update matrix Ã.
    The top n rows of Ã are zero (not stored explicitly).
    The top n rows of A_rem contain Householder reflectors.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update. Together with the top n rows of A_rem,
    this can be used to apply the block Householder transformation to other matrices.
    The number of rows depends on the block size and is architecture-dependent.
)doc");

    // Applying Householder transformations
    m.def(
        "update_apply_householder",
        [](c_matrix L, c_matrix A, c_matrix_cm W, c_matrix_cm B) {
            check_dim(L, A, W, B);
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::py::apply_householder(view(L̃), view(Ã), hyhound::Update{},
                                           view(W), view(B));
            return std::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "W"_a, "B"_a,
        R"doc(
Apply a block Householder transformation generated during a Cholesky
factorization update. Returns updated copies.

(L̃ Ã) = (L A) Q

where Q is the block Householder transformation represented by W and B.

Parameters
----------
L : (l × n), rectangular
    Matrix to apply the transformation to.

A : (l × m), rectangular
    Matrix to apply the transformation to.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update.

B : (k × m), rectangular
    The Householder reflector vectors generated during the Cholesky
    factorization update.

Returns
-------
L̃ : (l × n)
    The updated matrix L.

Ã : (l × m)
    The updated matrix A.
)doc");
    m.def(
        "downdate_apply_householder",
        [](c_matrix L, c_matrix A, c_matrix_cm W, c_matrix_cm B) {
            check_dim(L, A, W, B);
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::py::apply_householder(
                view(L̃), view(Ã), hyhound::Downdate{}, view(W), view(B));
            return std::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "W"_a, "B"_a,
        R"doc(
Apply a block Householder transformation generated during a Cholesky
factorization downdate. Returns updated copies.

(L̃ Ã) = (L A) Q

where Q is the block Householder transformation represented by W and B.

Parameters
----------
L : (l × n), rectangular
    Matrix to apply the transformation to.

A : (l × m), rectangular
    Matrix to apply the transformation to.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization downdate.

B : (k × m), rectangular
    The Householder reflector vectors generated during the Cholesky
    factorization downdate.

Returns
-------
L̃ : (l × n)
    The updated matrix L.

Ã : (l × m)
    The updated matrix A.
)doc");

    m.def(
        "update_apply_householder_sign",
        [](c_matrix L, c_matrix A, c_vector signs, c_matrix_cm W,
           c_matrix_cm B) {
            check_dim(L, A, W, B);
            if (A.shape(1) != signs.size())
                throw std::invalid_argument("len(signs) should be A.shape[1]");
            std::span signs_span{signs.data(), signs.shape(0)};
            if (!std::ranges::all_of(signs_span, [](T t) { return t == T{}; }))
                throw std::invalid_argument("signs should be +/- zero");
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::UpDowndate<T> sgn{signs_span};
            hyhound::py::apply_householder(view(L̃), view(Ã), sgn, view(W),
                                           view(B));
            return std::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "signs"_a, "W"_a, "B"_a,
        R"doc(
Apply a block Householder transformation generated during a Cholesky
factorization update with signed columns. Returns updated copies.

(L̃ Ã) = (L A) Q

where Q is the block Householder transformation represented by W, B and signs.

Parameters
----------
L : (l × n), rectangular
    Matrix to apply the transformation to.

A : (l × m), rectangular
    Matrix to apply the transformation to.

signs : m-vector
    Signs that determine whether a column of A was added (+0) or removed (-0).
    Values other than ±0 are not allowed.
    
W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update.

B : (k × m), rectangular
    The Householder reflector vectors generated during the Cholesky
    factorization update.

Returns
-------
L̃ : (l × n)
    The updated matrix L.

Ã : (l × m)
    The updated matrix A.
)doc");

    m.def(
        "update_apply_householder_diag",
        [](c_matrix L, c_matrix A, c_vector diag, c_matrix_cm W,
           c_matrix_cm B) {
            check_dim(L, A, W, B);
            if (A.shape(1) != diag.size())
                throw std::invalid_argument("len(diag) should be A.shape[1]");
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::py::apply_householder(view(L̃), view(Ã), d, view(W),
                                           view(B));
            return std::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "diag"_a, "W"_a, "B"_a,
        R"doc(
Apply a block Householder transformation generated during a Cholesky
factorization update with diagonal scaling. Returns updated copies.

(L̃ Ã) = (L A) Q

where Q is the block Householder transformation represented by W, B and diag.

Parameters
----------
L : (l × n), rectangular
    Matrix to apply the transformation to.

A : (l × m), rectangular
    Matrix to apply the transformation to.

diag : m-vector
    Scale factors corresponding to the columns of A used when generating the
    Householder transformation.

W : (r × n)
    The upper triangular Householder representations generated during the
    Cholesky factorization update.

B : (k × m), rectangular
    The Householder reflector vectors generated during the Cholesky
    factorization update.

Returns
-------
L̃ : (l × n)
    The updated matrix L.

Ã : (l × m)
    The updated matrix A.
)doc");
}

} // namespace

NB_MODULE(MODULE_NAME, m) {
    m.attr("__version__") = HYHOUND_VERSION_FULL;
    m.attr("build_time")  = HYHOUND_BUILD_TIME;
    m.attr("commit_hash") = HYHOUND_COMMIT_HASH;
#if HYHOUND_WITH_DOUBLE
    register_module<double>(m);
#endif
#if HYHOUND_WITH_FLOAT
    register_module<float>(m);
#endif
}
