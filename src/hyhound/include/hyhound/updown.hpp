#pragma once

#include <hyhound/config.hpp>
#include <span>
#include <type_traits>
#include <vector>

namespace hyhound {

/// Perform a factorization update, i.e. given the factorization of @f$ H @f$,
/// compute the factorization of @f$ H + A A^\top @f$.
struct Update {};
/// Perform a factorization downdate, i.e. given the factorization of @f$ H @f$,
/// compute the factorization of @f$ H - A A^\top @f$.
struct Downdate {};
/// Perform a factorization update or downdate, depending on the given signs,
/// i.e. given the factorization of @f$ H @f$, compute the factorization of
/// @f$ H + A S A^\top @f$, where @f$ S @f$ is a diagonal matrix with
/// @f$ S_{jj} = 1 @f$ if @p signs[j] is `+0.0`, and @f$ S_{jj} = -1 @f$ if
/// @p signs[j] is `-0.0`. Other values for @p signs are not allowed.
template <class T>
struct UpDowndate {
    std::span<const T> signs;
    [[nodiscard]] index_t size() const {
        return static_cast<index_t>(signs.size());
    }
};
/// Perform a factorization downdate or update, depending on the given signs,
/// i.e. given the factorization of @f$ H @f$, compute the factorization of
/// @f$ H - A S A^\top @f$, where @f$ S @f$ is a diagonal matrix with
/// @f$ S_{jj} = 1 @f$ if @p signs[j] is `+0.0`, and @f$ S_{jj} = -1 @f$ if
/// @p signs[j] is `-0.0`. Other values for @p signs are not allowed.
template <class T>
struct DownUpdate {
    std::span<const T> signs;
    [[nodiscard]] index_t size() const {
        return static_cast<index_t>(signs.size());
    }
};
/// Perform a factorization update or downdate with a general diagonal matrix,
/// i.e. given the factorization of @f$ H @f$, compute the factorization of
/// @f$ H + A D A^\top @f$.
template <class T>
struct DiagonalUpDowndate {
    std::span<const T> diag;
    [[nodiscard]] index_t size() const {
        return static_cast<index_t>(diag.size());
    }
};

template <class U>
UpDowndate(std::span<U>) -> UpDowndate<std::remove_const_t<U>>;
template <class U>
DownUpdate(std::span<U>) -> DownUpdate<std::remove_const_t<U>>;
template <class U>
DiagonalUpDowndate(std::span<U>) -> DiagonalUpDowndate<std::remove_const_t<U>>;

template <class U, class A>
UpDowndate(const std::vector<U, A> &) -> UpDowndate<U>;
template <class U, class A>
DownUpdate(const std::vector<U, A> &) -> DownUpdate<U>;
template <class U, class A>
DiagonalUpDowndate(const std::vector<U, A> &) -> DiagonalUpDowndate<U>;

} // namespace hyhound
