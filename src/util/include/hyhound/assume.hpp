#pragma once

#include <hyhound/config.hpp>
#include <guanaqo/assume.hpp>

/// @def HYHOUND_ASSUME(x)
/// Invokes undefined behavior if the expression @p x does not evaluate to true.
/// @throws std::logic_error in debug mode (when `NDEBUG` is not defined).

#if defined(NDEBUG) && !HYHOUND_VERIFY_ASSUMPTIONS
#define HYHOUND_ASSUME(x) GUANAQO_ASSUME(x)
#endif // defined(NDEBUG) && !HYHOUND_VERIFY_ASSUMPTIONS

#define HYHOUND_ASSERT(x) GUANAQO_ASSERT(x)

#ifndef HYHOUND_ASSUME
#define HYHOUND_ASSUME(x) HYHOUND_ASSERT(x)
#endif
