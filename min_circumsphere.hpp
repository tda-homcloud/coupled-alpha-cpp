#ifndef COUPLED_ALPHA_MIN_CIRCUMSPHERE
#define COUPLED_ALPHA_MIN_CIRCUMSPHERE

#include <Eigen/Dense>
#include "common.hpp"
#include "simplex.hpp"

namespace coupled_alpha {

template<int D>
std::pair<double, Eigen::Vector<double, D>>
MinimalCircumsphere(const Simplex& simplex, const std::vector<Eigen::Vector<double, D>>& coords) {
  using Vector = Eigen::Vector<double, D>;
  
  switch (simplex.Dim()) {
    case 0:
      return {0.0, coords[simplex[0]]};
    case 1:
      {
        auto& p0 = coords[simplex[0]];
        auto& p1 = coords[simplex[1]];
        return {(p0 - p1).norm() / 2.0, (p0 + p1) / 2.0};
      }
    case 2:
    case 3:
      {
        Eigen::MatrixXd A(simplex.Dim(), D);
        Eigen::VectorXd b(simplex.Dim());
        const Vector& p0 = coords[simplex[0]];

        for (int i = 1; i <= simplex.Dim(); ++i) {
          const Vector d = coords[simplex[i]] - p0;
          A.row(i - 1) = d.transpose();
          b[i - 1] = 0.5 * d.squaredNorm();
        }
        // std::cout << A << std::endl;
        // std::cout << b << std::endl;
        
        const Vector r = A.colPivHouseholderQr().solve(b);
        return {r.norm(), p0 + r};
      }
    default:
      assert(0);
  }
}

template<int D>
struct RelaxedFiltrationValue {
  Vectord<D> center;
  double r, rx, ry;

  static RelaxedFiltrationValue Compute(const Simplex& sx, const Simplex& sy,
                                        const std::vector<Vectord<D>>& coords) {
    if (sx.Empty()) {
      auto [r, c] = MinimalCircumsphere<D>(sy, coords);
      return {c, r, -1.0, r};
    }
    if (sy.Empty()) {
      auto [r, c] = MinimalCircumsphere<D>(sx, coords);
      return {c, r, r, -1.0};
    }

    Vectord<D> x0 = coords[sx[0]];
    Vectord<D> y0 = coords[sy[0]];
    int M = sx.Dim();
    int N = sy.Dim();
    
    if (N == 0 && M == 0) {
      double r = (x0 - y0).norm() / 2.0;
      return {(x0 + y0) / 2.0, r, r, r};
    }

    Eigen::MatrixXd Ap(N + M + 1, D);
    Eigen::VectorXd bp(N + M + 1);

    for (int m = 1; m <= M; ++m) {
      Ap.row(m - 1) = (coords[sx[m]] - x0).transpose();
      bp[m - 1] = 0.5 * (coords[sx[m]].squaredNorm() - x0.squaredNorm());
    }

    for (int n = 1; n <= N; ++n) {
      Ap.row(M + n - 1) = (coords[sy[n]] - y0).transpose();
      bp[M + n - 1] = 0.5 * (coords[sy[n]].squaredNorm() - y0.squaredNorm());
    }

    auto A = Ap.topLeftCorner(N + M, D);
    auto b = bp.topLeftCorner(N + M, 1);
    // std::cout << A << std::endl;
    // std::cout << b << std::endl;
    
    Vectord<D> dx = A.colPivHouseholderQr().solve(b - A * x0);
    Vectord<D> cx = x0 + dx;
    double rxx = dx.norm();
    double rxy = (cx - y0).norm();
    if (rxx >= rxy)
      return {cx, rxx, rxx, rxy};

    Vectord<D> dy = A.colPivHouseholderQr().solve(b - A * y0);
    Vectord<D> cy = y0 + dy;
    double ryy = dy.norm();
    double ryx = (cy - x0).norm();
    if (ryy >= ryx)
      return {cy, ryy, ryx, ryy};

    Ap.row(N + M) = (x0 - y0).transpose();
    bp[N + M] = 0.5 * (x0.squaredNorm() - y0.squaredNorm());
    Vectord<D> d = Ap.colPivHouseholderQr().solve(bp - Ap * x0);
    double r = d.norm();
    return {d + x0, r, r, r};
  }
};


} // namespace coupled_alpha

#endif
