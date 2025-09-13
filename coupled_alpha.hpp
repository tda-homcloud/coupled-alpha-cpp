#ifndef COUPLED_ALPHA_COUPLED_ALPHA
#define COUPLED_ALPHA_COUPLED_ALPHA

#include <vector>
#include <cstdint>

#include "common.hpp"
#include "simplex.hpp"
#include "min_circumsphere.hpp"

namespace coupled_alpha {

template<int D>
bool CoupledGabriel(
    const RelaxedFiltrationValue<D>& face_rfv,
    const Simplex& simplex,
    const std::vector<Vectord<D>>& coords,
    const std::vector<uint8_t>& levels) {

  auto [sx, sy] = simplex.Split(levels);
  for (int i = 0; i <= sx.Dim(); ++i) {
    if ((coords[sx[i]] - face_rfv.center).norm() <= face_rfv.rx - EPS)
      return false;
  }
  for (int i = 0; i <= sy.Dim(); ++i) {
    if ((coords[sy[i]] - face_rfv.center).norm() <= face_rfv.ry - EPS)
      return false;
  }
  return true;
}


} //namespace coupled_alpha


#endif
