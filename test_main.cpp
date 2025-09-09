#include <cmath>

#include "simplex.hpp"
#include "lifted_delaunay.hpp"
#include "min_circumsphere.hpp"

static void test_lifted_delaunay() {
  using V = Eigen::Vector2d;
  using Cell = coupled_alpha::LiftedDelaunay<2>::Cell;
  
  coupled_alpha::LiftedDelaunay<2> lf;
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

static void test_simplex() {
  using Simplex = coupled_alpha::Simplex;
  Simplex s = Simplex{0, 2, 4};
  std::cout << s << std::endl;
  std::cout << s.Dim() << std::endl;
}

static void test_min_circumsphere() {
  using Simplex = coupled_alpha::Simplex;
  using Vector2d = Eigen::Vector2d;
  using Vector3d = Eigen::Vector3d;
  using std::sqrt;
  
  std::vector<Vector2d> points2d = {
    Vector2d{0.0, 4.0}, Vector2d{2.0, 2.0},
    Vector2d{0.0, 0.0}, Vector2d{1.0, 0.0}, Vector2d{0.5, sqrt(3.0) / 2.0},
  };

  {
    auto [r, c] = coupled_alpha::MinimalCircumsphere<2>(Simplex{0}, points2d);
    assert(r == 0.0);
    assert((c - Vector2d{0.0, 4.0}).norm() < 1e-5);
  }

  {
    auto [r, c] = coupled_alpha::MinimalCircumsphere<2>(Simplex{0, 1}, points2d);
    assert(r == sqrt(2.0));
    assert((c - Vector2d{1.0, 3.0}).norm() < 1e-5);
  }
  {
    auto [r, c] = coupled_alpha::MinimalCircumsphere<2>(Simplex{2, 3, 4}, points2d);
    assert(std::fabs(r - 1 / sqrt(3.0)) < 1e-5);
    assert((c - Vector2d{0.5, 1.0 / (2 * sqrt(3.0))}).norm() < 1e-5);
  }

  std::vector<Vector3d> points3d = {
    Vector3d{0.0, 0.0, 0.0}, Vector3d{1.0, 0.0, 0.0}, Vector3d{0.5, sqrt(3.0) / 2.0, 0.0},
    Vector3d{0.5, sqrt(3.0) / 6.0, sqrt(6.0) / 3.0},
  };

  {
    auto [r, c] = coupled_alpha::MinimalCircumsphere<3>(Simplex{0, 1, 2, 3}, points3d);
    assert(std::fabs(r - sqrt(3.0 / 8.0)) < 1e-5);
    assert((c - Vector3d{0.5, 1.0 / (2 * sqrt(3.0)), sqrt(6.0) / 12}).norm() < 1e-5);
  }
}

int main(int argc, char** argv) {
  test_lifted_delaunay();
  test_simplex();
  test_min_circumsphere();
  return 0;
}
