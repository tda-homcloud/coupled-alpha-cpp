#include <CGAL/config.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <Eigen/Dense>

#include <cstdint>
#include <cassert>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <array>
#include <algorithm>
#include <initializer_list>

namespace coupled_alpha {

const double EPS  = 1e-10;

class Simplex {
 public:
  Simplex(): vertices_{UNUSED, UNUSED, UNUSED, UNUSED} {}

  Simplex(std::initializer_list<uint32_t> init) {
    std::copy(init.begin(), init.end(), vertices_);
    std::fill(vertices_ + init.size(), vertices_ + 4, UNUSED);
  }

  std::size_t dim() const {
    for (std::size_t i = 0; i < 4; ++i) {
      if (vertices_[i] == UNUSED)
        return i - 1;
    }
    return 3;
  }

  uint32_t operator[](size_t k) const {
    assert(k < 4);
    assert(vertices_[k] != UNUSED);
    return vertices_[k];
  }
  
  size_t hash() const {
    return static_cast<size_t>(vertices_[0] ^ vertices_[1] ^ vertices_[2] ^ vertices_[3]);
  }
 private:
  static const uint32_t UNUSED = 0xffffffff;
  uint32_t vertices_[4];
};

class SimplexHash {
  size_t operator()(const Simplex& simplex) { return simplex.hash(); }
};

template<int D>
class CoupledAlphaComplex {
  static_assert(D == 2 || D == 3);
 public:
  using K = CGAL::Epick_d<CGAL::Dimension_tag<D + 1>>;
  using T = CGAL::Delaunay_triangulation<K>;
  using Point = typename T::Point;
  using Vertex_handle = typename T::Vertex_handle;

  using Cell = std::array<uint32_t, D + 2>;
  using Vector = Eigen::Vector<double, D>;

  CoupledAlphaComplex(): dt_(D + 1), n_(0) {}
  
  void InsertCoord(const Point& point) {
    if constexpr (D == 2) {
      coords_.emplace_back(point[0], point[1]);
    } else {
      coords_.emplace_back(point[0], point[1], point[2]);
    }
  }
  
  void Insert(const Point& point) {
    assert(point[D] == 0.0 || point[D] == 1.0);
    InsertCoord(point);
    labels_.push_back(static_cast<uint8_t>(point[D]));
    
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

  std::pair<double, Vector>
  MinimalCircumsphere(const Simplex& simplex) {
    switch (simplex.dim()) {
      case 0:
        return std::make_pair(0.0, coords_[simplex[0]]);
      case 1: {
        const Vector& x = coords_[simplex[0]];
        const Vector& y = coords_[simplex[1]];
        return std::make_pair((x - y).norm() / 2, (x + y) / 2);
      }
      case 2:
      case 3: {
        Eigen::MatrixXd A(simplex.dim(), D);
        Vector b;
        
        const Vector& p0 = coords_[simplex[0]];
        for (int i = 1; i <= simplex.dim(); ++i) {
          const Vector d = p0 - coords_[simplex[i]];
          A.row(i - 1) = d.transpose();
          b[i - 1] = 0.5 * d.squaredNorm();
        }

        Vector r = A.colPivHouseholderQr().solve(b);
        return {p0 + r, r.norm()};
      }
      default:
        assert(0);
    }
  }

  static void PrintCoord(const Vector& c) {
    std::cout << "[";
    for (size_t i = 0; i < D; ++i)
      std::cout << c[i] << ", ";
    std::cout << "]" << std::endl;
  }

  static void PrintCell(const Cell& cell) {
    std::cout << "{";
    for (std::size_t i = 0; i < D + 2; ++i) {
      std::cout << cell[i];
      if (i != D + 1) {
        std::cout << ", ";
      }
    }
    std::cout << "}" << std::endl;
  }

  void Print() {
    std::cout << "Points" << std::endl;
    for (const auto& c: coords_) {
      PrintCoord(c);
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
  std::vector<Vector> coords_;
  std::vector<uint8_t> labels_;
  T dt_;
  uint32_t n_;
  Vertex_handle last_handle_;
  std::unordered_map<Vertex_handle, uint32_t> vertex_ids_;
  std::vector<Cell> cells_;

 public:
  class RelaxedFiltrationValue {
   public:
    Point center_;
    double r_, r_0, r_1;
  };
};

  
} // namespace coupled_alpha

