#ifndef COUPLED_ALPHA_SIMPLEX_HPP
#define COUPLED_ALPHA_SIMPLEX_HPP

#include <cstdint>
#include <cassert>
#include <initializer_list>
#include <array>
#include <unordered_set>
#include <ostream>
#include <utility>
#include <functional>
#include <iterator>
#include <iostream>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/container/static_vector.hpp>

#include "common.hpp"

namespace coupled_alpha {

class SimplexFaceIterator;

class Simplex {  
 public:
  using Container = boost::container::static_vector<uint32_t, 5>;
  
  friend class SimplexFaceIterator;
  
  Simplex() {}

  Simplex(std::initializer_list<uint32_t> init): vertices_(init) {}

  inline bool operator==(const Simplex& other) const {
    return vertices_ == other.vertices_;
  }
  
  template<typename Iterator>
  Simplex(Iterator begin, Iterator end): vertices_(begin, end) {}
  
  inline void Append(uint32_t n) {
    vertices_.push_back(n);
  }

  inline bool Empty() const { return vertices_.empty(); }

  inline const Container& Vertices() const { return vertices_; }

  inline uint16_t Dim() const {
    assert(vertices_.size() > 0);
    return vertices_.size() - 1;
  }

  inline uint32_t operator[](std::size_t i) const {
    return vertices_[i];
  }

  inline std::pair<Simplex, Simplex> Split(const std::vector<uint8_t>& labels) const {
    Simplex s0, s1;
    for (size_t i = 0; i < vertices_.size(); ++i) {
      if (labels[vertices_[i]] == 0) {
        s0.Append(vertices_[i]);
      } else {
        s1.Append(vertices_[i]);
      }
    }
    return {s0, s1};
  }

  inline std::size_t hash() const {
    std::size_t hash_value = 0;
    for (uint32_t value: vertices_)
      hash_value ^= value;
    return hash_value;
  }

  inline std::vector<Simplex> ProperFaces() const {
    if (Dim() == 0)
      return std::vector<Simplex>();
    
    std::vector<Simplex> faces;
    Simplex face(vertices_.begin() + 1, vertices_.end());
    
    for (size_t i = 0; i < vertices_.size() - 1; ++i) {
      faces.push_back(face);
      face.vertices_[i] = vertices_[i];
    }
    faces.push_back(face);

    return faces;
  }
  
  inline friend std::ostream& operator<<(std::ostream& os, const Simplex& simplex) {
    if (simplex.vertices_.empty()) {
      os << "Simplex{}";
    } else {
      os << "Simplex{";
      for (size_t i = 0; i < simplex.vertices_.size() - 1; ++i) {
        os << simplex.vertices_[i] << ", ";
      }
      os << simplex.vertices_[simplex.vertices_.size() - 1] << "}";
    }
    return os;
  }

 private:
  Container vertices_;

};

struct SimplexHash {
  typedef std::size_t result_type;
  std::size_t operator()(const Simplex& simplex) const { return simplex.hash(); }
};


template<int D>
class CellsToSimpices {
 public:
  using Simplices = std::array<std::unordered_set<Simplex, SimplexHash>, D + 2>;

  static Simplices Compute(const std::vector<Cell<D>>& cells) {
    Simplices simplices;
    CellsToSimpices builder(simplices);
    
    for (const auto& cell: cells) {
      Simplex top_simplex(cell.begin(), cell.end());
      builder.visit(top_simplex);
    }
    
    return simplices;
  }

 private:
  Simplices& simplices_;
  CellsToSimpices(Simplices& simplices): simplices_(simplices) {}

  void visit(const Simplex& simplex) {
    if (simplices_[simplex.Dim()].count(simplex) == 0) {
      simplices_[simplex.Dim()].insert(simplex);
      for (const auto& face: simplex.ProperFaces())
        visit(face);
    }
  }
};

} // namespace coupled_alpha

#endif
