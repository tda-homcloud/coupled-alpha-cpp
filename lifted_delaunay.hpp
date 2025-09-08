#ifndef COUPLED_ALPHA_LIFTED_DELAUNAY_HPP
#define COUPLED_ALPHA_LIFTED_DELAUNAY_HPP

#include <array>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <Eigen/Dense>

namespace coupled_alpha {

template<int D>
class LiftedDelaunay {
  static_assert(D == 2 || D == 3);

 public:
  using K = CGAL::Epick_d<CGAL::Dimension_tag<D + 1>>;
  using T = CGAL::Delaunay_triangulation<K>;
  using Vertex_handle = typename T::Vertex_handle;
  using Vector = Eigen::Vector<double, D>;
  using LPoint = typename T::Point;
  using Cell = std::array<uint32_t, D + 2>;
  
 private:
  T dt_;
  uint32_t n_;
  Vertex_handle last_handle_;
  std::unordered_map<Vertex_handle, uint32_t> vertex_ids_;
  
 public:
  LiftedDelaunay(): dt_(D + 1), n_(0) {}

  void Insert(const Vector& coord, uint8_t level) {
    LPoint lifted_point = Lift(coord, level);
    
    if (last_handle_ == Vertex_handle()) {
      last_handle_ = dt_.insert(lifted_point);
    } else {
      last_handle_ = dt_.insert(lifted_point, last_handle_);
    }
    vertex_ids_[last_handle_] = n_;
    ++n_;
  }

  LPoint Lift(const Vector& c, uint8_t level) {
    if constexpr (D == 2) {
      return LPoint(c[0], c[1], level);
    } else {
      return LPoint(c[0], c[1], c[2], level);
    }
  }

  std::vector<Cell> Cells() {
    std::vector<Cell> cells;

    for (auto cit = dt_.full_cells_begin(); cit != dt_.full_cells_end(); ++cit) {
      if (dt_.is_infinite(cit))
        continue;

      Cell& cell = cells.emplace_back();
      int k = 0;

      for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit, ++k) {
        cell[k] = vertex_ids_[*vit];
      }
      std::sort(cell.begin(), cell.end());
    }

    return cells;
  }
};

} // namespace coupled_alpha

#endif
