#pragma once

#include <hyhound/config.hpp>
#include <guanaqo/lut.hpp>
#include <type_traits>

namespace hyhound {

template <index_t I>
using index_constant = std::integral_constant<index_t, I>;

/// Returns a 2D array of the form:
///
/// ~~~
/// {{   f(0, 0),     f(0, 1),    ...,   f(0, C - 1)   },
///  {   f(1, 0),     f(1, 1),    ...,   f(1, C - 1)   },
///  {     ...,         ...,      ...,       ...       },
///  { f(R - 1, 0), f(R - 1, 1)}, ..., f(R - 1, C - 1) }}
/// ~~~
///
/// The argument @p f should be a function (or callable) that accepts two
/// arguments of type @ref index_constant.
template <int R, int C, class F>
consteval auto make_2d_lut(F f) {
    return guanaqo::make_2d_lut<index_t, R, C>(std::forward<F>(f));
}

/// Returns an array of the form:
///
/// ~~~
/// {   f(0),     f(1),    ...,   f(C - 1)   }
/// ~~~
///
/// The argument @p f should be a function (or callable) that accepts an
/// argument of type @ref index_constant.
template <int N, class F>
consteval auto make_1d_lut(F f) {
    return guanaqo::make_1d_lut<index_t, N>(std::forward<F>(f));
}

} // namespace hyhound
