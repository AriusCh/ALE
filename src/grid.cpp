#include "grid.hpp"

#include <cmath>

Grid::Grid(GridType type_, size_t sizeX_, size_t sizeY_)
    : sizeX(sizeX_), sizeY(sizeY_), type(type_) {}

GridALE::GridALE(const std::vector<std::vector<double>> &x,
                 const std::vector<std::vector<double>> &y,
                 const std::vector<std::vector<double>> &rho,
                 const std::vector<std::vector<double>> &p,
                 const std::vector<std::vector<double>> &u,
                 const std::vector<std::vector<double>> &v,
                 std::unique_ptr<EOS> &&eos)
    : Grid(GridType::eGridALE, rho.size(), rho[0].size()),
      x(x),
      y(y),
      rho(rho),
      p(p),
      u(u),
      v(v),
      eos(std::move(eos)) {}
