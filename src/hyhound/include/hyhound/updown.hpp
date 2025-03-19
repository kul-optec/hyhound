#pragma once

#include <hyhound/config.hpp>
#include <span>

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
struct UpDowndate {
    std::span<const real_t> signs;
    [[nodiscard]] index_t size() const {
        return static_cast<index_t>(signs.size());
    }
};
/// Perform a factorization downdate or update, depending on the given signs,
/// i.e. given the factorization of @f$ H @f$, compute the factorization of
/// @f$ H - A S A^\top @f$, where @f$ S @f$ is a diagonal matrix with
/// @f$ S_{jj} = 1 @f$ if @p signs[j] is `+0.0`, and @f$ S_{jj} = -1 @f$ if
/// @p signs[j] is `-0.0`. Other values for @p signs are not allowed.
struct DownUpdate {
    std::span<const real_t> signs;
    [[nodiscard]] index_t size() const {
        return static_cast<index_t>(signs.size());
    }
};
/// Perform a factorization update or downdate with a general diagonal matrix,
/// i.e. given the factorization of @f$ H @f$, compute the factorization of
/// @f$ H + A D A^\top @f$.
struct DiagonalUpDowndate {
    std::span<const real_t> diag;
    [[nodiscard]] index_t size() const {
        return static_cast<index_t>(diag.size());
    }
};

} // namespace hyhound
