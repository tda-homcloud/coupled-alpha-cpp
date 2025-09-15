
#include <cstdint>
#include <cstdlib>

#include "simplex.hpp"
#include "lifted_delaunay.hpp"

using namespace coupled_alpha;

int main(int argc, char** argv) {
  using V = Eigen::Vector2d;
  using Cell = LiftedDelaunay<2>::Cell;
  
  LiftedDelaunay<2> lf;
  lf.insert(V{0.0, 0.00}, 0);
  lf.insert(V{1.0, 0.00}, 0);
  lf.insert(V{0.5, 0.72}, 0);
  lf.insert(V{0.0, 0.72}, 1);
  lf.insert(V{1.0, 0.72}, 1);
  lf.insert(V{0.5, 0.00}, 1);

  std::vector<Cell> cells = lf.cells();
  for (const auto& cell: cells)
    std::cout << Simplex(cell.begin(), cell.end()) << std::endl;

  return 0;
}
