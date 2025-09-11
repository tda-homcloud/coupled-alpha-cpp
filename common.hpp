#ifndef COUPLED_ALPHA_COMMON
#define COUPLED_ALPHA_COMMON

#include <Eigen/Dense>

namespace coupled_alpha {

static const double EPS = 1e-10;

template<int D>
using Vectord = Eigen::Vector<double, D>;

} // namespace coupled_alpha

#endif
