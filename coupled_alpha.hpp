#ifndef COUPLED_ALPHA_COUPLED_ALPHA
#define COUPLED_ALPHA_COUPLED_ALPHA

#include <vector>
#include <cstdint>

#include "common.hpp"
#include "simplex.hpp"
#include "min_circumsphere.hpp"
#include "lifted_delaunay.hpp"

namespace coupled_alpha {

template<int D>
class CoupledAlpha {
 public:
  CoupledAlpha(
      const std::vector<Vectord<D>>& coords,
      const std::vector<uint8_t>& levels):
      coords_(coords), levels_(levels) {}
  
  void compute() {
    LiftedDelaunay<D> ldelaunay;
    for (size_t i = 0; i < coords_.size(); ++i)
      ldelaunay.insert(coords_[i], levels_[i]);

    const auto cells = ldelaunay.cells();
    const auto simplices = CellsToSimpices<D>::compute(cells);
    std::array<std::unordered_map<Simplex, RelaxedFiltrationValue<D>, SimplexHash>, D + 2> rfvs;
    
    for (int d = 0; d <= D + 1; ++d) {
      for (const auto& simplex: simplices[d]) {
        auto [sx, sy] = simplex.split(levels_);
        rfvs[d][simplex] = RelaxedFiltrationValue<D>::compute(sx, sy, coords_);
        std::cout << simplex.join(" ") << std::setprecision(7) << rfvs[d][simplex].r << std::endl;
      }
    }
  }
  
  // bool is_coupled_gabriel(const Simplex& face, const Simplex& simplex) const {
  // }
 private:
  const std::vector<Vectord<D>>& coords_;
  const std::vector<uint8_t>& levels_;
};
  
// template<int D>
// bool is_coupled_gabriel(
//     const RelaxedFiltrationValue<D>& face_rfv,
//     const Simplex& simplex,
//     const std::vector<Vectord<D>>& coords,
//     const std::vector<uint8_t>& levels) {

//   auto [sx, sy] = simplex.Split(levels);
//   for (auto v: sx.Vertices()) {
//     if ((coords[v] - face_rfv.center).norm() <= face_rfv.rx - EPS)
//       return false;
//   }
//   for (auto v: sy.Vertices()) {
//     if ((coords[v] - face_rfv.center).norm() <= face_rfv.ry - EPS)
//       return false;
//   }
//   return true;
// }

// template<int D>
// std::array<std::unordered_map<Simplex, double, SimplexHash>, D + 2>
// CoupledAlpha(const std::vector<Vectord<D>>& coords,
//              const std::vector<uint8_t>& levels){
//   LiftedDelaunay<D> ldelaunay;
//   for (size_t i = 0; i < coords.size(); ++i)
//     ldelaunay.Insert(coords[i], levels[i]);
  
//   auto cells = ldelaunay.Cells();
//   auto simplices = CellsToSimpices<D>::Compute(cells);
//   std::array<std::unordered_map<Simplex, RelaxedFiltrationValue<D>, SimplexHash>, D + 2> rfvs;
//   std::array<std::unordered_map<Simplex, double, SimplexHash>, D + 2> values;
  
//   for (int k = 1; k <= D + 1; ++k) {
//     for (const auto& P: simplices[k]) {
//       auto [sx, sy] = P.Split(levels);
//       rfvs[k][P] = RelaxedFiltrationValue<D>::Compute(sx, sy, coords);
//     }
//   }
  
//   for (int k = D + 1; k >= 1; --k) {
//     for (const auto& P: simplices[k]) {
//       if (values[k].count(P) == 0) {
//         values[k].insert(std::make_pair(P, rfvs[k][P].r));
//       }
//       for (const auto& Q: P.ProperFaces()) {
//         if (values[k - 1].count(Q) == 1) {
//           values[k - 1][Q] = std::min(values[k - 1][Q], values[k][P]);
//         } else {
//           if (!CoupledGabriel(rfvs[k - 1][Q], P, coords, levels)) {
//             values[k - 1].insert(std::make_pair(Q, values[k][P]));
//           }
//         }
//       }
//     }
//   }
//   for (const auto& pt: simplices[0]) {
//     values[0][pt] = 0;
//   }
//   return values;
// }

  
} //namespace coupled_alpha


#endif
