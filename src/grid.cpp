#include "grid.hpp"

#include <cmath>

Grid::Grid(GridType type_, size_t sizeX_, size_t sizeY_)
    : sizeX(sizeX_), sizeY(sizeY_), type(type_) {}

GridALE::GridALE(size_t sizeX_, size_t sizeY_,
                 const std::unique_ptr<Polygon> &polygon)

    : Grid(GridType::eGridALE, sizeX_, sizeY_) {
  polygon->generateMesh(x, y);

  auto bcs = rect->getBoundaries();

  bottomBcs = std::make_unique<Boundary>(bcs[0]);
  rightBcs = std::make_unique<Boundary>(bcs[1]);
  topBcs = std::make_unique<Boundary>(bcs[2]);
  leftBcs = std::make_unique<Boundary>(bcs[3]);
}
