#include "coupled_alpha.hpp"

using Point = typename coupled_alpha::CoupledAlphaComplex<2>::Point;

int main(int argc, char** argv) {

  coupled_alpha::CoupledAlphaComplex<2> cac;
  
  cac.Insert(Point{0.0, 0.0, 0.0});
  cac.Insert(Point{1.0, 0.0, 0.0});
  cac.Insert(Point{0.5, 0.72, 0.0});
  cac.Insert(Point{0.0, 0.72, 1.0});
  cac.Insert(Point{1.0, 0.72, 1.0});
  cac.Insert(Point{0.5, 0.0, 1.0});
  cac.ComputeDelaunay();
  cac.Print();

  return 0;
}
