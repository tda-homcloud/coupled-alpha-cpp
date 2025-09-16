#include "coupled_alpha.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace coupled_alpha;

static std::vector<std::string> split_line(const std::string& str) {
  std::vector<std::string> ret;

  std::string token;
  std::stringstream ss(str);
  while (ss >> token)
    ret.push_back(std::move(token));

  return ret;
}

template<int D>
std::pair<std::vector<uint8_t>, std::vector<Vectord<D>>>
load_txt(std::ifstream& file) {
  std::string line;
  std::vector<Vectord<D>> coords;
  std::vector<uint8_t> levels;
  
  while (std::getline(file, line)) {
    if (line[0] == '#')
      continue;
    auto tokens = split_line(line);
    assert(D == static_cast<int>(tokens.size()) - 1);

    levels.push_back(std::stoi(tokens[0]));
    Vectord<D> coord;

    for (size_t i = 0; i < D; ++i)
      coord[i] = std::stod(tokens[i + 1]);
    
    coords.push_back(std::move(coord));
  }
  return {levels, coords};
}

template<int D>
void run(std::ifstream& f) {
  auto [levels, coords] = load_txt<D>(f);
  CoupledAlpha<D> coupled_alpha(coords, levels);
  const auto filtration_values = coupled_alpha.compute();

  for (size_t d = 0; d <= D + 1; ++d) {
    for (const auto& [simplex, value]: filtration_values[d]) {
      std::cout << simplex.join(" ") << " "
                << std::fixed << std::setprecision(7) << value << std::endl;
    }
  }
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cout << "Usage: ./coupled_alpha DIM FNAME" << std::endl;
    std::exit(1);
  }
  int dim = std::stoi(std::string(argv[1]));
  std::ifstream f(argv[2]);

  if (dim == 2) {
    run<2>(f);
  } else if (dim == 3) {
    run<3>(f);
  } else {
    std::cerr << "Unsupported dimension: " << dim << std::endl;
  }

  return 0;
}
