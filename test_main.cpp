#include <cmath>
#include <vector>
#include <utility>
#include <cstdint>
#include <cstdlib>

#include "simplex.hpp"
#include "lifted_delaunay.hpp"
#include "min_circumsphere.hpp"

using namespace coupled_alpha;

static void TestLiftedDelaunay() {
  using V = Eigen::Vector2d;
  using Cell = LiftedDelaunay<2>::Cell;
  
  LiftedDelaunay<2> lf;
  lf.Insert(V{0.0, 0.00}, 0);
  lf.Insert(V{1.0, 0.00}, 0);
  lf.Insert(V{0.5, 0.72}, 0);
  lf.Insert(V{0.0, 0.72}, 1);
  lf.Insert(V{1.0, 0.72}, 1);
  lf.Insert(V{0.5, 0.00}, 1);

  std::vector<Cell> cells = lf.Cells();

  for (const Cell& cell: cells) {
    std::cout << "{"
              << cell[0] << ", "
              << cell[1] << ", "
              << cell[2] << ", "
              << cell[3] << "}" << std::endl;
  }
}

static void TestSimplex() {
  using Simplex = Simplex;
  std::vector<uint8_t> labels = {0, 1, 1, 0, 0, 0};
  Simplex s = Simplex{0, 2, 4};
  
  std::cout << s << std::endl;
  std::cout << s.Dim() << std::endl;
  auto [s0, s1] = s.Split(labels);
  std::cout << s0 << " " << s1 << std::endl;
}

static void TestMinimalCircumsphere() {
  using Simplex = Simplex;
  using Vector2d = Eigen::Vector2d;
  using Vector3d = Eigen::Vector3d;
  using std::sqrt;
  
  std::vector<Vector2d> points2d = {
    Vector2d{0.0, 4.0}, Vector2d{2.0, 2.0},
    Vector2d{0.0, 0.0}, Vector2d{1.0, 0.0}, Vector2d{0.5, sqrt(3.0) / 2.0},
  };

  {
    auto [r, c] = MinimalCircumsphere<2>(Simplex{0}, points2d);
    assert(r == 0.0);
    assert((c - Vector2d{0.0, 4.0}).norm() < 1e-5);
  }

  {
    auto [r, c] = MinimalCircumsphere<2>(Simplex{0, 1}, points2d);
    assert(r == sqrt(2.0));
    assert((c - Vector2d{1.0, 3.0}).norm() < 1e-5);
  }
  {
    auto [r, c] = MinimalCircumsphere<2>(Simplex{2, 3, 4}, points2d);
    assert(std::fabs(r - 1 / sqrt(3.0)) < 1e-5);
    assert((c - Vector2d{0.5, 1.0 / (2 * sqrt(3.0))}).norm() < 1e-5);
  }

  std::vector<Vector3d> points3d = {
    Vector3d{0.0, 0.0, 0.0}, Vector3d{1.0, 0.0, 0.0}, Vector3d{0.5, sqrt(3.0) / 2.0, 0.0},
    Vector3d{0.5, sqrt(3.0) / 6.0, sqrt(6.0) / 3.0},
  };

  {
    auto [r, c] = MinimalCircumsphere<3>(Simplex{0, 1, 2, 3}, points3d);
    assert(std::fabs(r - sqrt(3.0 / 8.0)) < 1e-5);
    assert((c - Vector3d{0.5, 1.0 / (2 * sqrt(3.0)), sqrt(6.0) / 12}).norm() < 1e-5);
  }
}

void AssertAlmostEq(double expected, double actual, double eps=1e-5) {
  if (fabs(expected - actual) >= eps) {
    std::cout << "AssertAlmostEq Expected: " << expected << ", Actual: " << actual << std::endl;
    std::exit(1);
  }
}

template<int D>
void AssertRelaxedFiltrationValue(
    const Eigen::Vector<double, D>& expected_center,
    double expected_r, double expected_rx, double expected_ry,
    const RelaxedFiltrationValue<D>& rfv) {
  if ((expected_center - rfv.center).norm() > 1e-5) {
    std::cout << "Expected: " << expected_center.transpose()
              << " Actual: " << rfv.center.transpose() << std::endl;
    std::exit(1);
  }
  AssertAlmostEq(expected_r, rfv.r);
  AssertAlmostEq(expected_rx, rfv.rx);
  AssertAlmostEq(expected_ry, rfv.ry);
}

static void TestRelaxedFiltrationValue() {
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  using std::sqrt;
  
  {
    std::vector<Vector2d> points2d = {
      Vector2d{0, 4}, Vector2d{2, 2},
      Vector2d{0, 0}, Vector2d{1, 0}, Vector2d{2, 1}, Vector2d{2, 5}, Vector2d{2, 3},
      Vector2d{2, 0}, Vector2d{5, 1}, Vector2d{5, 5}, Vector2d{3, 1}, Vector2d{3, 5},
    };
    
    RelaxedFiltrationValue<2> rfv;

    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{0, 1}, Simplex{}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{1.0, 3.0}, sqrt(2.0), sqrt(2.0), -1, rfv);
    
    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{}, Simplex{0, 1}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{1.0, 3.0}, sqrt(2.0), -1.0, sqrt(2.0), rfv);

    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{0}, Simplex{1}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{1.0, 3.0}, sqrt(2.0), sqrt(2.0), sqrt(2.0), rfv);
    
    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{2, 3}, Simplex{4}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{0.5, 1.0}, 1.5, sqrt(1.25), 1.5, rfv);

    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{4}, Simplex{2, 3}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{0.5, 1.0}, 1.5, 1.5, sqrt(1.25), rfv);
    
    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{2}, Simplex{4, 6}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{0.25, 2.0}, sqrt(4.0625), sqrt(4.0625), sqrt(4.0625), rfv);

    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{2, 7}, Simplex{8, 9}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{1, 3}, sqrt(20), sqrt(10), sqrt(20), rfv);

    rfv = RelaxedFiltrationValue<2>::Compute(Simplex{2, 7}, Simplex{10, 11}, points2d);
    AssertRelaxedFiltrationValue(Vector2d{1, 3}, sqrt(10), sqrt(10), sqrt(8), rfv);
  }

  {
    std::vector<Vector3d> points3d = {
      Vector3d{0, 0, 0}, Vector3d{2, 0, 0},
      Vector3d{5, 1, 2}, Vector3d{5, 5, 2},
      Vector3d{3, 1, 6}, Vector3d{3, 3, 6},
    };

    RelaxedFiltrationValue<3> rfv;
    
    rfv = RelaxedFiltrationValue<3>::Compute(Simplex{0, 1}, Simplex{2, 3}, points3d);
    AssertRelaxedFiltrationValue(Vector3d{1, 3, 2}, sqrt(20), sqrt(14), sqrt(20), rfv);

    rfv = RelaxedFiltrationValue<3>::Compute(Simplex{0, 1}, Simplex{4, 5}, points3d);
    AssertRelaxedFiltrationValue(Vector3d{1, 2, 3}, sqrt(14), sqrt(14), sqrt(14), rfv);

  }
}

                                  

int main(int argc, char** argv) {
  TestLiftedDelaunay();
  TestSimplex();
  TestMinimalCircumsphere();
  TestRelaxedFiltrationValue();
  return 0;
}
