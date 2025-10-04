#include <hyhound/householder-updowndate.hpp>
#include <hyhound/updown.hpp>
#include <hyhound-version.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <memory>

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
            hyhound::update_cholesky(view(L), view(A), hyhound::Update{});
        },
        "L"_a.noconvert(), "A"_a.noconvert());
    m.def(
        "downdate_cholesky_inplace",
        [](matrix_cm L, matrix_cm A) {
            hyhound::update_cholesky(view(L), view(A), hyhound::Downdate{});
        },
        "L"_a.noconvert(), "A"_a.noconvert());
    m.def(
        "update_cholesky_sign_inplace",
        [](matrix_cm L, matrix_cm A, c_vector signs) {
            hyhound::UpDowndate<T> sgn{std::span{signs.data(), signs.shape(0)}};
            hyhound::update_cholesky(view(L), view(A), sgn);
        },
        "L"_a.noconvert(), "A"_a.noconvert(), "signs"_a);
    m.def(
        "update_cholesky_diag_inplace",
        [](matrix_cm L, matrix_cm A, c_vector diag) {
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::update_cholesky(view(L), view(A), d);
        },
        "L"_a.noconvert(), "A"_a.noconvert(), "diag"_a);
    // Returning copies
    m.def(
        "update_cholesky",
        [](c_matrix L, c_matrix A) {
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::update_cholesky(view(L̃), view(Ã), hyhound::Update{});
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a);
    m.def(
        "downdate_cholesky",
        [](c_matrix L, c_matrix A) {
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::update_cholesky(view(L̃), view(Ã), hyhound::Downdate{});
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a);
    m.def(
        "update_cholesky_sign",
        [](c_matrix L, c_matrix A, c_vector signs) {
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::UpDowndate<T> sgn{std::span{signs.data(), signs.shape(0)}};
            hyhound::update_cholesky(view(L̃), view(Ã), sgn);
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "signs"_a);
    m.def(
        "update_cholesky_diag",
        [](c_matrix L, c_matrix A, c_vector diag) {
            auto L̃ = copy(L), Ã = copy(A);
            hyhound::DiagonalUpDowndate<T> d{
                std::span{diag.data(), diag.shape(0)}};
            hyhound::update_cholesky(view(L̃), view(Ã), d);
            return nb::make_tuple(std::move(L̃), std::move(Ã));
        },
        "L"_a, "A"_a, "diag"_a);
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
