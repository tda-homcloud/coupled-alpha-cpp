#include <cmath>
#include <cassert>

#include "simplex.hpp"

using namespace coupled_alpha;

int main(int argc, char** argv) {
  std::vector<uint32_t> vec = {0, 2, 4};
  Simplex s(vec.begin(), vec.end());

  std::cout << s << std::endl;
  std::cout << Simplex{0, 2, 4} << std::endl;
  assert((Simplex{0, 2, 4} == s));
  auto faces = s.faces();
  assert((faces.size() == 3));

  auto [s0, s1] = s.split(std::vector<uint8_t>{0, 0, 0, 1, 1, 1});
  assert((s0 == Simplex{0, 2}));
  assert((s1 == Simplex{4}));

  std::vector<Cell<3>> cells = {Cell<3>{0, 1, 2, 3, 4}, Cell<3>{0, 1, 2, 3, 5}};
  const auto simplices = CellsToSimpices<3>::compute(cells);
  assert(simplices[4].size() == 2);
  assert(simplices[3].size() == 9);
  assert(simplices[2].size() == 4 + 12);
  assert(simplices[1].size() == 6 + 4 * 2);
  assert(simplices[0].size() == 6);
  
  return 0;
}
