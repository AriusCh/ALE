#include "grid.hpp"

#include <cmath>

Grid::Grid(GridType type_, size_t sizeX_, size_t sizeY_)
    : sizeX(sizeX_), sizeY(sizeY_), type(type_) {}

GridALE::GridALE(size_t sizeX_, size_t sizeY_,
                 std::unique_ptr<PrimitiveRectangle> rect)
    : Grid(GridType::eGridALE, sizeX_, sizeY_) {
  std::array<double, 4> xv;
  std::array<double, 4> yv;
  rect->getVertices(xv, yv);

  x.reserve(sizeX + 1);
  y.reserve(sizeX + 1);
  double dxi = (xv[3] - xv[0]) / sizeX;
  double dyi = (yv[3] - yv[0]) / sizeY;
  double dxj = (xv[1] - xv[0]) / sizeX;
  double dyj = (yv[1] - yv[0]) / sizeY;
  for (int i = 0; i < sizeX + 1; i++) {
    std::vector<double> tmpX;
    std::vector<double> tmpY;
    tmpX.reserve(sizeY + 1);
    tmpY.reserve(sizeY + 1);
    for (int j = 0; j < sizeY + 1; j++) {
      tmpX.push_back(xv[0] + dxi * i + dxj * j);
      tmpY.push_back(yv[0] + dyi * i + dyj * j);
    }
    x.push_back(tmpX);
    y.push_back(tmpY);
  }

  auto bcs = rect->getBoundaries();

  bottomBcs = std::make_unique<Boundary>(bcs[0]);
  rightBcs = std::make_unique<Boundary>(bcs[1]);
  topBcs = std::make_unique<Boundary>(bcs[2]);
  leftBcs = std::make_unique<Boundary>(bcs[3]);
}
