#include <cassert>

#include "common.hpp"
#include "min_circumsphere.hpp"

using namespace coupled_alpha;

static void assert_almost_eq(double expected, double actual, double eps=1e-5) {
  if (fabs(expected - actual) >= eps) {
    std::cout << "Expected: " << expected << ", Actual: " << actual << std::endl;
    std::exit(1);
  }
}

template<int D>
void assert_relaxed_filtration_value(
    const Vectord<D>& expected_center,
    double expected_r, double expected_rx, double expected_ry,
    const RelaxedFiltrationValue<D>& rfv) {
  if ((expected_center - rfv.center).norm() > 1e-5) {
    std::cout << "Expected: " << expected_center.transpose()
              << " Actual: " << rfv.center.transpose() << std::endl;
    std::exit(1);
  }
  assert_almost_eq(expected_r, rfv.r);
  assert_almost_eq(expected_rx, rfv.rx);
  assert_almost_eq(expected_ry, rfv.ry);
}

int main(int argc, char** argv) {
  using Eigen::Vector2d;
  
  std::vector<Vector2d> points2d = {
    Vector2d{0, 4}, Vector2d{2, 2},
    Vector2d{0, 0}, Vector2d{1, 0}, Vector2d{2, 1}, Vector2d{2, 5}, Vector2d{2, 3},
    Vector2d{2, 0}, Vector2d{5, 1}, Vector2d{5, 5}, Vector2d{3, 1}, Vector2d{3, 5},
  };
    
  RelaxedFiltrationValue<2> rfv;

  rfv = RelaxedFiltrationValue<2>::compute(Simplex{0, 1}, Simplex{}, points2d);
  assert_relaxed_filtration_value(Vector2d{1.0, 3.0}, 2, 2, -1, rfv);
    
  rfv = RelaxedFiltrationValue<2>::compute(Simplex{}, Simplex{0, 1}, points2d);
  assert_relaxed_filtration_value(Vector2d{1.0, 3.0}, 2, -1.0, 2, rfv);

  rfv = RelaxedFiltrationValue<2>::compute(Simplex{0}, Simplex{1}, points2d);
  assert_relaxed_filtration_value(Vector2d{1.0, 3.0}, 2, 2, 2, rfv);
    
  rfv = RelaxedFiltrationValue<2>::compute(Simplex{2, 3}, Simplex{4}, points2d);
  assert_relaxed_filtration_value(Vector2d{0.5, 1.0}, 2.25, 1.25, 2.25, rfv);

  rfv = RelaxedFiltrationValue<2>::compute(Simplex{4}, Simplex{2, 3}, points2d);
  assert_relaxed_filtration_value(Vector2d{0.5, 1.0}, 2.25, 2.25, 1.25, rfv);
    
  rfv = RelaxedFiltrationValue<2>::compute(Simplex{2}, Simplex{4, 6}, points2d);
  assert_relaxed_filtration_value(Vector2d{0.25, 2.0}, 4.0625, 4.0625, 4.0625, rfv);

  rfv = RelaxedFiltrationValue<2>::compute(Simplex{2, 7}, Simplex{8, 9}, points2d);
  assert_relaxed_filtration_value(Vector2d{1, 3}, 20, 10, 20, rfv);

  rfv = RelaxedFiltrationValue<2>::compute(Simplex{2, 7}, Simplex{10, 11}, points2d);
  assert_relaxed_filtration_value(Vector2d{1, 3}, 10, 10, 8, rfv);

  return 0;
}
