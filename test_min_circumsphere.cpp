#include <cmath>

#include "min_circumsphere.hpp"

using namespace coupled_alpha;

template<int D>
void assert_min_circumsphere(
    double expected_radius,
    const Vectord<D>& expected_center,
    const Simplex& simplex,
    const std::vector<Vectord<D>>& coords) {
  auto [r, c] = min_circumsphere(simplex, coords);
  assert(std::fabs(r - expected_radius) < 1e-8);
  assert((c - expected_center).norm() < 1e-8);
}
    
int main(int argc, char** argv) {
  using Vector2d = Eigen::Vector2d;
  using Vector3d = Eigen::Vector3d;
  using std::sqrt;
  
  {
    std::vector<Vector2d> points2d = {
      Vector2d{0.0, 4.0}, Vector2d{2.0, 2.0},
      Vector2d{0.0, 0.0}, Vector2d{1.0, 0.0}, Vector2d{0.5, sqrt(3.0) / 2.0},
    };
    
    assert_min_circumsphere(0.0, Vector2d{0, 4}, Simplex{0}, points2d);
    assert_min_circumsphere(sqrt(2), Vector2d{1, 3}, Simplex{0, 1}, points2d);
    assert_min_circumsphere(1 / sqrt(3), Vector2d{0.5, 1 / (2 * sqrt(3))}, Simplex{2, 3, 4}, points2d);
  }
  {
    std::vector<Vector3d> points3d = {
      Vector3d{0, 0, 0}, Vector3d{1, 0, 0}, Vector3d{0.5, sqrt(3) / 2, 0}, 
      Vector3d{0.5, sqrt(3) / 6, sqrt(6) / 3},
    };
    assert_min_circumsphere(sqrt(3.0 / 8),
                            Vector3d{0.5, 1 / (2 * sqrt(3)), sqrt(6) / 12},
                            Simplex{0, 1, 2, 3}, points3d);
  }
}
