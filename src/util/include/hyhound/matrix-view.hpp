#pragma once

#include <hyhound/config.hpp>
#include <guanaqo/mat-view.hpp>

namespace hyhound {

using RealMatrixView        = guanaqo::MatrixView<const real_t, index_t>;
using MutableRealMatrixView = guanaqo::MatrixView<real_t, index_t>;

} // namespace hyhound
