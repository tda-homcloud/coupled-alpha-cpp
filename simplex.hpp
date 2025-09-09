#ifndef COUPLED_ALPHA_SIMPLEX_HPP
#define COUPLED_ALPHA_SIMPLEX_HPP

#include <cstdint>
#include <cassert>
#include <initializer_list>
#include <array>
#include <ostream>

namespace coupled_alpha {

class Simplex {
 private:
  uint16_t num_verteices_;
  std::array<uint32_t, 4> vertices_;
  
 public:
  Simplex(): num_verteices_(0) {}

  Simplex(std::initializer_list<uint32_t> init) {
    assert(init.size() <= 4);
    num_verteices_ = init.size();
    std::copy(init.begin(), init.end(), vertices_.begin());
  }

  void Append(uint32_t n) {
    assert(num_verteices_ < 4);
    vertices_[num_verteices_] = n;
    ++num_verteices_;
  }
  
  inline uint16_t Dim() const {
    assert(num_verteices_ > 0);
    return num_verteices_ - 1;
  }
  
  uint32_t operator[](std::size_t i) const {
    assert(i < num_verteices_);
    return vertices_[i];
  }
  
  std::size_t hash() const {
    std::size_t hash_value = 0;
    for (uint16_t i = 0; i < num_verteices_; ++i)
      hash_value ^= vertices_[i];
    return hash_value;
  }

  friend std::ostream& operator<<(std::ostream& os, const Simplex& simplex) {
    os << "Simplex{";
    for (uint16_t i = 0; i < simplex.num_verteices_ - 1; ++i) {
      os << simplex.vertices_[i] << ", ";
    }
    os << simplex.vertices_[simplex.num_verteices_ - 1] << "}";
    return os;
  }
};

class SimplexHash {
  std::size_t operator()(const Simplex& simplex) { return simplex.hash(); }
};

} // namespace coupled_alpha

#endif
