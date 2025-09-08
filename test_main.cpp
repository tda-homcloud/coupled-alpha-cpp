#include "lifted_delaunay.hpp"

void test_lifted_delaunay() {
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


int main(int argc, char** argv) {
  test_lifted_delaunay();
  return 0;
}
