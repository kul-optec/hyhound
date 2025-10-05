#include <hyhound/householder-updowndate.hpp>
#include <hyhound/updown.hpp>
#include <hyhound-version.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <memory>
#include <stdexcept>

namespace nb = nanobind;
using namespace nb::literals;

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

template <class T>
void register_module(nb::module_ &m) {
    using c_matrix = nb::ndarray<const T, nb::ndim<2>, nb::device::cpu>;
    using c_vector =
        nb::ndarray<const T, nb::ndim<1>, nb::device::cpu, nb::any_contig>;
    using matrix_cm =
        nb::ndarray<T, nb::ndim<2>, nb::device::cpu, nb::f_contig>;
    // In-place
    m.def(
        "update_cholesky_inplace",
        [](matrix_cm L, matrix_cm A) {
            check_dim(L, A);
            hyhound::update_cholesky(view(L), view(A), hyhound::Update{});
        },
        "L"_a.noconvert(), "A"_a.noconvert(),
        "Cholesky factorization update. Overwrites its arguments.\n\n"
        "L̃L̃ᵀ + ÃÃᵀ = LLᵀ + AAᵀ\n\n"
        ":param L: k × n matrix, lower-trapezoidal, Fortran order. On entry, "
        "the original Cholesky factor L. "
        "On exit, it contains the updated Cholesky factor L̃.\n"
        ":param A: k × m matrix, rectangular, Fortran order. On entry, the "
        "update matrix A. "
        "On exit, it contains the k-n bottom rows of the remaining update "
        "matrix Ã (the top n rows of Ã are implicitly zero). The top n rows of "
        "this matrix are overwritten by the Householder reflectors used during "
        "the update, and are generally not useful.\n");
    m.def(
        "downdate_cholesky_inplace",
        [](matrix_cm L, matrix_cm A) {
            check_dim(L, A);
            hyhound::update_cholesky(view(L), view(A), hyhound::Downdate{});
        },
        "L"_a.noconvert(), "A"_a.noconvert(),
        "Cholesky factorization downdate. Overwrites its arguments.\n\n"
        "L̃L̃ᵀ - ÃÃᵀ = LLᵀ - AAᵀ\n\n"
        ":param L: k × n matrix, lower-trapezoidal, Fortran order. On entry, "
        "the original Cholesky factor L. "
        "On exit, it contains the updated Cholesky factor L̃.\n"
        ":param A: k × m matrix, rectangular, Fortran order. On entry, the "
        "downdate matrix A. "
        "On exit, it contains the k-n bottom rows of the remaining downdate "
        "matrix Ã (the top n rows of Ã are implicitly zero). The top n rows of "
        "this matrix are overwritten by the Householder reflectors used during "
        "the downdate, and are generally not useful.\n");
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
            hyhound::update_cholesky(view(L), view(A), sgn);
        },
        "L"_a.noconvert(), "A"_a.noconvert(), "signs"_a,
        "Cholesky factorization update. Overwrites its arguments.\n\n"
        "L̃L̃ᵀ + ÃSÃᵀ = LLᵀ + ASAᵀ, "
        "with ``S = np.diag(np.copysign(np.ones(m), signs)``, and where "
        "``signs`` contains ±0.\n\n"
        ":param L: k × n matrix, lower-trapezoidal, Fortran order. On entry, "
        "the original Cholesky factor L. "
        "On exit, it contains the updated Cholesky factor L̃.\n"
        ":param A: k × m matrix, rectangular, Fortran order. On entry, the "
        "update matrix A. "
        "On exit, it contains the k-n bottom rows of the remaining update "
        "matrix Ã (the top n rows of Ã are implicitly zero). The top n rows of "
        "this matrix are overwritten by the Householder reflectors used during "
        "the update, and are generally not useful.\n"
        ":param signs: m-vector. Signs that determine whether a column of A is "
        "added (+0) or removed (-0). Values other than ±0 are not allowed.");
    m.def(
        "update_cholesky_diag_inplace",
        [](matrix_cm L, matrix_cm A, c_vector diag) {
            check_dim(L, A);
            if (A.shape(1) != diag.size())
                throw std::invalid_argument("len(diag) should be A.shape[1]");
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::update_cholesky(view(L), view(A), d);
        },
        "L"_a.noconvert(), "A"_a.noconvert(), "diag"_a,
        "Cholesky factorization update. Overwrites its arguments.\n\n"
        "L̃L̃ᵀ + ÃDÃᵀ = LLᵀ + ADAᵀ, "
        "with ``D = np.diag(diag)``.\n\n"
        ":param L: k × n matrix, lower-trapezoidal, Fortran order. On entry, "
        "the original Cholesky factor L. "
        "On exit, it contains the updated Cholesky factor L̃.\n"
        ":param A: k × m matrix, rectangular, Fortran order. On entry, the "
        "update matrix A. "
        "On exit, it contains the k-n bottom rows of the remaining update "
        "matrix Ã (the top n rows of Ã are implicitly zero). The top n rows of "
        "this matrix are overwritten by the Householder reflectors used during "
        "the update, and are generally not useful.\n"
        ":param diag: m-vector. Scale factors corresponding to the columns of "
        "A.\n");
    // Returning copies
    m.def(
        "update_cholesky",
        [](c_matrix L, c_matrix A) {
            check_dim(L, A);
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::update_cholesky(view(L̃), view(Ã), hyhound::Update{});
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a,
        "Cholesky factorization update. Returns updated copies.\n\n"
        "L̃L̃ᵀ + ÃÃᵀ = LLᵀ + AAᵀ\n\n"
        ":param L: k × n matrix, lower-trapezoidal. The original Cholesky "
        "factor.\n"
        ":param A: k × m matrix, rectangular. The update matrix.\n"
        ":return: Tuple (L̃, Ã). L̃ is the updated Cholesky factor. Ã contains "
        "the k-n bottom rows of the remaining update matrix. "
        "The top n rows of Ã are overwritten by Householder reflectors and are "
        "generally not useful.\n");
    m.def(
        "downdate_cholesky",
        [](c_matrix L, c_matrix A) {
            check_dim(L, A);
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::update_cholesky(view(L̃), view(Ã), hyhound::Downdate{});
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a,
        "Cholesky factorization downdate. Returns updated copies.\n\n"
        "L̃L̃ᵀ - ÃÃᵀ = LLᵀ - AAᵀ\n\n"
        ":param L: k × n matrix, lower-trapezoidal. The original Cholesky "
        "factor.\n"
        ":param A: k × m matrix, rectangular. The downdate matrix.\n"
        ":return: Tuple (L̃, Ã). L̃ is the updated Cholesky factor. Ã contains "
        "the k-n bottom rows of the remaining update matrix. "
        "The top n rows of Ã are overwritten by Householder reflectors and are "
        "generally not useful.\n");
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
            hyhound::UpDowndate<T> sgn{signs_span};
            hyhound::update_cholesky(view(L̃), view(Ã), sgn);
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "signs"_a,
        "Cholesky factorization update. Returns updated copies.\n\n"
        "L̃L̃ᵀ + ÃSÃᵀ = LLᵀ + ASAᵀ, "
        "with ``S = np.diag(np.copysign(np.ones(m), signs))``, and where "
        "``signs`` contains ±0.\n\n"
        ":param L: k × n matrix, lower-trapezoidal. The original Cholesky "
        "factor.\n"
        ":param A: k × m matrix, rectangular. The update matrix.\n"
        ":param signs: m-vector. Signs that determine whether a column of A is "
        "added (+0) or removed (-0). "
        "Values other than ±0 are not allowed.\n"
        ":return: Tuple (L̃, Ã). L̃ is the updated Cholesky factor. Ã contains "
        "the k-n bottom rows of the remaining update matrix. "
        "The top n rows of Ã are overwritten by Householder reflectors and are "
        "generally not useful.\n");
    m.def(
        "update_cholesky_diag",
        [](c_matrix L, c_matrix A, c_vector diag) {
            check_dim(L, A);
            if (A.shape(1) != diag.size())
                throw std::invalid_argument("len(diag) should be A.shape[1]");
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::update_cholesky(view(L̃), view(Ã), d);
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "diag"_a,
        "Cholesky factorization update. Returns updated copies.\n\n"
        "L̃L̃ᵀ + ÃDÃᵀ = LLᵀ + ADAᵀ, "
        "with ``D = np.diag(diag)``.\n\n"
        ":param L: k × n matrix, lower-trapezoidal. The original Cholesky "
        "factor.\n"
        ":param A: k × m matrix, rectangular. The update matrix.\n"
        ":param diag: m-vector. Scale factors corresponding to the columns of "
        "A.\n"
        ":return: Tuple (L̃, Ã). L̃ is the updated Cholesky factor. Ã contains "
        "the k-n bottom rows of the remaining update matrix. "
        "The top n rows of Ã are overwritten by Householder reflectors and are "
        "generally not useful.\n");
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
