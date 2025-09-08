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

#include "simplex.hpp"
#include "lifted_delaunay.hpp"

namespace coupled_alpha {


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

