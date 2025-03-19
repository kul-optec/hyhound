#pragma once

#include <hyhound/config.hpp>
#include <experimental/simd>
#include <cstddef>

namespace hyhound::micro_kernels {

namespace stdx = std::experimental;

template <class T, class OuterStrideT = index_t>
struct mat_access_impl {
    using value_type = T;
    value_type *data;
    [[no_unique_address]] OuterStrideT outer_stride;
    static constexpr ptrdiff_t inner_stride = 1;
    static constexpr bool transpose         = false;

    [[gnu::always_inline]] value_type &operator()(index_t r,
                                                  index_t c) const noexcept {
        ptrdiff_t i0 = transpose ? c : r;
        ptrdiff_t i1 = transpose ? r : c;
        return data[i0 * inner_stride +
                    i1 * static_cast<ptrdiff_t>(outer_stride)];
    }
    template <class Simd>
    [[gnu::always_inline]] Simd load(index_t r, index_t c) const noexcept {
        return Simd{&operator()(r, c), stdx::element_aligned};
    }
    template <class Simd>
    [[gnu::always_inline]] void store(Simd x, index_t r,
                                      index_t c) const noexcept
        requires(!std::is_const_v<T>)
    {
        x.copy_to(&operator()(r, c), stdx::element_aligned);
    }
    template <class Simd, class Align>
    [[gnu::always_inline]] Simd load(index_t r, index_t c,
                                     Align align) const noexcept {
        return Simd{&operator()(r, c), align};
    }
    template <class Simd, class Align>
    [[gnu::always_inline]] void store(Simd x, index_t r, index_t c,
                                      Align align) const noexcept
        requires(!std::is_const_v<T>)
    {
        x.copy_to(&operator()(r, c), align);
    }
    template <class Self>
    [[gnu::always_inline]] constexpr Self block(this const Self &self,
                                                index_t r, index_t c) noexcept {
        return {&self(r, c), self.outer_stride};
    }
    template <class Self>
    [[gnu::always_inline]] constexpr Self middle_rows(this const Self &self,
                                                      index_t r) noexcept {
        return {&self(r, 0), self.outer_stride};
    }
    template <class Self>
    [[gnu::always_inline]] constexpr Self middle_cols(this const Self &self,
                                                      index_t c) noexcept {
        return {&self(0, c), self.outer_stride};
    }

    [[gnu::always_inline]] constexpr mat_access_impl(
        value_type *data = nullptr, OuterStrideT outer_stride = {}) noexcept
        : data{data}, outer_stride{outer_stride} {}
    [[gnu::always_inline]] constexpr mat_access_impl(
        const guanaqo::MatrixView<T, index_t> &o) noexcept
        : data{o.data},
          outer_stride{o.outer_stride * static_cast<index_t>(inner_stride)} {}
    [[gnu::always_inline]] constexpr mat_access_impl(
        const guanaqo::MatrixView<std::remove_const_t<T>, index_t> &o) noexcept
        requires(std::is_const_v<T>)
        : data{o.data},
          outer_stride{o.outer_stride * static_cast<index_t>(inner_stride)} {}
    [[gnu::always_inline]] constexpr mat_access_impl(
        const mat_access_impl<std::remove_const_t<T>> &o) noexcept
        requires(std::is_const_v<T>)
        : data{o.data}, outer_stride{o.outer_stride} {}
    [[gnu::always_inline]] constexpr mat_access_impl(const mat_access_impl &o) =
        default;
};

struct matrix_accessor : mat_access_impl<const real_t> {
    using mat_access_impl<const real_t>::mat_access_impl;
};

struct mut_matrix_accessor : mat_access_impl<real_t> {
    using mat_access_impl<real_t>::mat_access_impl;
};

} // namespace hyhound::micro_kernels
