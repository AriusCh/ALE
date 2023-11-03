#include "grid.hpp"

#include <cmath>

GridALE::GridALE(int sizeX_, int sizeY_, double xmin, double xmax, double ymin,
                 double ymax, std::function<double(double, double)> uInit,
                 std::function<double(double, double)> vInit,
                 std::function<double(double, double)> rhoInit,
                 std::function<double(double, double)> pInit,
                 std::shared_ptr<EOS> eos)
    : sizeX(sizeX_),
      sizeY(sizeY_),
      x(sizeX_ + 1, std::vector<double>(sizeY_ + 1)),
      y(sizeX_ + 1, std::vector<double>(sizeY_ + 1)),
      u(sizeX_ + 1, std::vector<double>(sizeY_ + 1)),
      v(sizeX_ + 1, std::vector<double>(sizeY_ + 1)),
      m(sizeX_ + 1, std::vector<double>(sizeY_ + 1)),
      rho(sizeX_, std::vector<double>(sizeY_)),
      p(sizeX_, std::vector<double>(sizeY_)),
      E(sizeX_, std::vector<double>(sizeY_)),
      eos(std::move(eos)) {
  // Init coordinates
  double dx = (xmax - xmin) / sizeX_;
  double dy = (ymax - ymin) / sizeY_;
  for (int i = 0; i < sizeX_ + 1; i++) {
    for (int j = 0; j < sizeY_ + 1; j++) {
      x[i][j] = xmin + dx * i;
      y[i][j] = ymin + dy * j;
    }
  }

  // Init velocities
  for (int i = 0; i < sizeX_ + 1; i++) {
    for (int j = 0; j < sizeY_ + 1; j++) {
      u[i][j] = uInit(x[i][j], y[i][j]);
      v[i][j] = vInit(x[i][j], y[i][j]);
    }
  }

  // Init densities and pressures
  for (int i = 0; i < sizeX_; i++) {
    for (int j = 0; j < sizeY_; j++) {
      double x_ =
          0.25 * (x[i + 1][j] + x[i + 1][j + 1] + x[i][j + 1] + x[i][j]);
      double y_ =
          0.25 * (y[i + 1][j] + y[i + 1][j + 1] + y[i][j + 1] + y[i][j]);
      rho[i][j] = rhoInit(x_, y_);
      p[i][j] = pInit(x_, y_);
    }
  }

  // Init auxulary values
  populateCellEnergy();
  populateNodeMass();
}
double GridALE::getV(int i, int j, std::vector<std::vector<double>> &x,
                     std::vector<std::vector<double>> &y) {
  double xij = x[i][j], xipj = x[i + 1][j], xipjp = x[i + 1][j + 1],
         xijp = x[i][j + 1];
  double yij = y[i][j], yipj = y[i + 1][j], yipjp = y[i + 1][j + 1],
         yijp = y[i][j + 1];
  return 0.5 * ((xij - xipjp) * (yipj - yijp) + (xipj - xijp) * (yipjp - yij));
}
void GridALE::populateCellEnergy() {
  // Populate energy
  for (int i = 0; i < sizeX; i++) {
    for (int j = 0; j < sizeY; j++) {
      double U = u[i + 1][j] * u[i + 1][j] + u[i + 1][j + 1] * u[i + 1][j + 1] +
                 u[i][j + 1] * u[i][j + 1] + u[i][j];
      double V = v[i + 1][j] * v[i + 1][j] + v[i + 1][j + 1] * v[i + 1][j + 1] +
                 v[i][j + 1] * v[i][j + 1] + v[i][j];
      E[i][j] = 0.125 * (U + V) + eos->gete(rho[i][j], p[i][j]);
    }
  }
}
void GridALE::populateNodeMass() {
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
  for (int i = 0; i < sizeX + 1; i++) {
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
      m[i][j] = 0.25 * (Mimjm + Mijm + Mij + Mimj);
    }
  }
}
