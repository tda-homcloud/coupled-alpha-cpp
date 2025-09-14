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

template<int D>
std::array<std::unordered_map<Simplex, double, SimplexHash>, D + 2>
CoupledAlpha(const std::vector<Vectord<D>>& coords,
             const std::vector<uint8_t>& levels){
  LiftedDelaunay<D> ldelaunay;
  for (size_t i = 0; i < coords.size(); ++i)
    ldelaunay.Insert(coords[i], levels[i]);
  
  auto cells = ldelaunay.Cells();
  auto simplices = CellsToSimpices<D>::Compute(cells);
  std::array<std::unordered_map<Simplex, RelaxedFiltrationValue<D>, SimplexHash>, D + 2> rfvs;
  std::array<std::unordered_map<Simplex, double, SimplexHash>, D + 2> values;
  
  for (int k = 1; k <= D + 1; ++k) {
    for (const auto& P: simplices[k]) {
      auto [sx, sy] = P.Split(levels);
      rfvs[k][P] = RelaxedFiltrationValue<D>::Compute(sx, sy, coords);
    }
  }
  
  for (int k = D + 1; k >= 1; --k) {
    for (const auto& P: simplices[k]) {
      if (values[k].count(P) == 0) {
        values[k].insert(std::make_pair(P, rfvs[k][P].r));
      }
      for (const auto& Q: P.ProperFaces()) {
        if (values[k - 1].count(Q) == 1) {
          values[k - 1][Q] = std::min(values[k - 1][Q], values[k][P]);
        } else {
          if (!CoupledGabriel(rfvs[k - 1][Q], P, coords, levels)) {
            values[k - 1][Q] = values[k][P];
          }
        }
      }
    }
  }
  return values;
}

  
} //namespace coupled_alpha


#endif
