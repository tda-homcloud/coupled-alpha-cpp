#ifndef COUPLED_ALPHA_SIMPLEX_HPP
#define COUPLED_ALPHA_SIMPLEX_HPP

#include <cstdint>
#include <initializer_list>

namespace coupled_alpha {

class Simplex {
 private:
  uint32_t vertices_[4];
  
  static const uint32_t UNUSED = 0xffffffff;
  
 public:
  Simplex(): vertices_{UNUSED, UNUSED, UNUSED, UNUSED} {}

  Simplex(std::initializer_list<uint32_t> init) {
    assert(init.size() <= 4);
    std::copy(init.begin(), init.end(), vertices_);
    std::fill(vertices_ + init.size(), vertices_ + 4, UNUSED);
  }

  std::size_t dim() const {
    for (std::size_t i = 1; i < 4; ++i) {
      if (vertices_[i] == UNUSED)
        return i - 1;
    }
    return 3;
  }

  uint32_t operator[](size_t i) const {
    assert(i < 4);
    assert(vertices_[i] != UNUSED);
    return vertices_[i];
  }
  
  size_t hash() const {
    return static_cast<size_t>(vertices_[0] ^ vertices_[1] ^ vertices_[2] ^ vertices_[3]);
  }
};

class SimplexHash {
  size_t operator()(const Simplex& simplex) { return simplex.hash(); }
};

} // namespace coupled_alpha

#endif
