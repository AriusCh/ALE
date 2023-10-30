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
      uNext(u),
      vNext(v),
      eos(std::move(eos)) {
  PopulateCellEnergy();
  PopulateNodeMass();
}
std::unique_ptr<Boundary> GridALE::createExternalBoundary(
    BoundaryType type, ExternalBoundarySide side) {
  return std::make_unique<ExternalBoundary>(type, side, x, y, rho, p, m, uNext,
                                            vNext);
}
double GridALE::getV(int i, int j, std::vector<std::vector<double>> &x,
                     std::vector<std::vector<double>> &y) const {
  double xij = x[i][j], xipj = x[i + 1][j], xipjp = x[i + 1][j + 1],
         xijp = x[i][j + 1];
  double yij = y[i][j], yipj = y[i + 1][j], yipjp = y[i + 1][j + 1],
         yijp = y[i][j + 1];
  return 0.5 * ((xij - xipjp) * (yipj - yijp) + (xipj - xijp) * (yipjp - yij));
}
void GridALE::PopulateCellEnergy() {
  E.clear();
  // Populate energy
  E.reserve(sizeX);
  for (int i = 0; i < sizeX; i++) {
    std::vector<double> ETemp;
    ETemp.reserve(sizeY);
    for (int j = 0; j < sizeY; j++) {
      double U = u[i + 1][j] * u[i + 1][j] + u[i + 1][j + 1] * u[i + 1][j + 1] +
                 u[i][j + 1] * u[i][j + 1] + u[i][j];
      double V = v[i + 1][j] * v[i + 1][j] + v[i + 1][j + 1] * v[i + 1][j + 1] +
                 v[i][j + 1] * v[i][j + 1] + v[i][j];
      ETemp.push_back(0.125 * (U + V) + eos->gete(rho[i][j], p[i][j]));
    }
    E.push_back(ETemp);
  }
}
void GridALE::PopulateNodeMass() {
  m.clear();
  // Temporary vector for cell mass
  std::vector<std::vector<double>> M;
  M.reserve(sizeX);
  for (int i = 0; i < sizeX; i++) {
    std::vector<double> MTemp;
    MTemp.reserve(sizeY);
    for (int j = 0; j < sizeY; j++) {
      MTemp.push_back(rho[i][j] * getV(i, j, x, y));
    }
    M.push_back(MTemp);
  }
  // Populate node mass
  m.reserve(sizeX + 1);
  for (int i = 0; i < sizeX + 1; i++) {
    std::vector<double> mTemp;
    mTemp.reserve(sizeY + 1);
    for (int j = 0; j < sizeY + 1; j++) {
      double Mimjm = 0.;
      double Mijm = 0.;
      double Mij = 0.;
      double Mimj = 0.;
      if (i != 0 && j != 0) {
        Mimjm = M[i - 1][j - 1];
      }
      if (i != sizeX && j != 0) {
        Mijm = M[i][j - 1];
      }
      if (i != sizeX && j != sizeY) {
        Mij = M[i][j];
      }
      if (i != 0 && j != sizeY) {
        Mimj = M[i - 1][j];
      }
      mTemp.push_back(0.25 * (Mimjm + Mijm + Mij + Mimj));
    }
    m.push_back(mTemp);
  }
}
