#ifndef COUPLED_ALPHA_MIN_CIRCUMSPHERE
#define COUPLED_ALPHA_MIN_CIRCUMSPHERE

#include <Eigen/Dense>

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
        
        Vector r = A.colPivHouseholderQr().solve(b);
        return {r.norm(), p0 + r};
      }
    default:
      assert(0);
  }
}
                    
} // namespace coupled_alpha

#endif
