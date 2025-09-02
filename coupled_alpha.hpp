#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include <cstdint>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <array>
#include <algorithm>


namespace coupled_alpha {

template<int dim>
class CoupledAlphaComplex {
 public:
  using K = CGAL::Epick_d<CGAL::Dimension_tag<dim + 1>>;
  using T = CGAL::Delaunay_triangulation<K>;
  using Point = typename T::Point;
  using Vertex_handle = typename T::Vertex_handle;

  using Cell = std::array<uint32_t, dim + 2>;
  
  CoupledAlphaComplex(): dt_(dim + 1), n_(0) {}

  void Insert(const Point& point) {
    assert(point[dim] == 0.0 || point[dim] == 1.0);
    points_.push_back(point);
    labels_.push_back(static_cast<uint8_t>(point[dim]));
    
    if (Vertex_handle() == last_handle_) {
      last_handle_ = dt_.insert(point);
    } else {
      last_handle_ = dt_.insert(point, last_handle_);
    }
    vertex_ids_[last_handle_] = n_;
    ++n_;
  }

  void ComputeDelaunay() {
    for (auto cit = dt_.full_cells_begin(); cit != dt_.full_cells_end(); ++cit) {
      if (dt_.is_infinite(cit))
        continue;
      auto& cell = cells_.emplace_back();

      int k = 0;
      for (auto vit = cit->vertices_begin(); vit != cit->vertices_end(); ++vit, ++k) {
        cell[k] = vertex_ids_[*vit];
      }
      std::sort(cell.begin(), cell.end());
    }
  }

  static void PrintPoint(const Point& point) {
    std::cout << "[";
    for (int i = 0; i < dim; ++i) {
      std::cout << point[i];
      if (i != dim - 1)
        std::cout << ", ";
    }
    std::cout << "](" << static_cast<int>(point[dim]) << ")" << std::endl;
  }
  
  static void PrintCell(const Cell& cell) {
    std::cout << "{";
    for (std::size_t i = 0; i < dim + 2; ++i) {
      std::cout << cell[i];
      if (i != dim + 1) {
        std::cout << ", ";
      }
    }
    std::cout << "}" << std::endl;
  }

  void Print() {
    std::cout << "Points" << std::endl;
    for (const auto& point: points_) {
      PrintPoint(point);
    }
    std::cout << "Labels" << std::endl;
    for (uint8_t label: labels_) {
      std::cout << static_cast<int>(label) << std::endl;
    }
    std::cout << "Cells" << std::endl;
    for (const auto& cell: cells_) {
      PrintCell(cell);
    }
  }
  
 private:
  std::vector<Point> points_;
  std::vector<uint8_t> labels_;
  T dt_;
  uint32_t n_;
  Vertex_handle last_handle_;
  std::unordered_map<Vertex_handle, uint32_t> vertex_ids_;
  std::vector<Cell> cells_;
};

  
} // namespace coupled_alpha

