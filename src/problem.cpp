#include "problem.hpp"

#include <cassert>

Problem::Problem(ProblemType type_, AxisymmetryType symType_)
    : type(type_), symType(symType_) {}

ProblemType Problem::getType() const { return type; }
AxisymmetryType Problem::getSymmetryType() const { return symType; }

std::unique_ptr<Grid> Problem::createGrid(size_t nx, size_t ny) const {
  const auto &polygons = calcRegion->getPolygons();

  for (size_t i = 0; i < polygons.size(); i++) {
    const auto &polygon = polygons[i];
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> v;
    polygon->generateMesh(x, y);
    initializeConditions[i](x, y, rho, p, u, v);
  }
}

RiemannProblem1D::RiemannProblem1D(double xmin, double xmax, double tmax,
                                   double rhoL, double uL, double pL,
                                   double rhoR, double uR, double pR,
                                   double spl, double gamma)
    : Problem(ProblemType::eRiemannProblem1D,
              AxisymmetryType::eNonSymmetrical) {
  assert(spl >= xmin && spl <= xmax);

  createGeometry(xmin, xmax);
  createBoundaryTypes();
  createInitialConditions(rhoL, pL, uL, rhoR, pR, uR, spl);
  createEOSes(gamma);
}
void RiemannProblem1D::createGeometry(double xmin, double xmax) {
  std::unique_ptr<Polygon> rect =
      std::make_unique<Rectangle>(xmin, 0.0, xmax - xmin, 1.0);
  calcRegion = std::make_unique<Region2D>(std::move(rect));
}
void RiemannProblem1D::createBoundaryTypes() {
  leftBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
  topBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
  rightBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
  bottomBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
}
void RiemannProblem1D::createInitialConditions(double rhoL, double pL,
                                               double uL, double rhoR,
                                               double pR, double uR,
                                               double spl) {
  auto initialCond = [rhoL, pL, uL, rhoR, pR, uR, spl](
                         const std::vector<std::vector<double>> &x,
                         const std::vector<std::vector<double>> &y,
                         std::vector<std::vector<double>> &rho,
                         std::vector<std::vector<double>> &p,
                         std::vector<std::vector<double>> &u,
                         std::vector<std::vector<double>> &v) {
    rho.clear();
    p.clear();
    u.clear();
    v.clear();
    rho.reserve(x.size() - 1);
    p.reserve(x.size() - 1);
    for (int i = 0; i < x.size() - 1; i++) {
      std::vector<double> rhoTmp;
      std::vector<double> pTmp;
      rhoTmp.reserve(y.size() - 1);
      pTmp.reserve(y.size() - 1);
      for (int j = 0; j < y.size() - 1; j++) {
        double x_ = 0.25 * (x[i + 1][j] + x[i][j + 1] + x[i][j + 1] + x[i][j]);
        if (x_ <= spl) {
          rhoTmp.push_back(rhoL);
          pTmp.push_back(pL);
        } else {
          rhoTmp.push_back(rhoR);
          pTmp.push_back(pR);
        }
      }
      rho.push_back(rhoTmp);
      p.push_back(pTmp);
    }

    u.reserve(x.size());
    v.reserve(y.size());
    for (int i = 0; i < x.size(); i++) {
      std::vector<double> uTmv;
      std::vector<double> vTmv;
      uTmv.reserve(y.size());
      vTmv.reserve(y.size());
      for (int j = 0; j < y.size(); j++) {
        double x_ = 0.25 * (x[i + 1][j] + x[i][j + 1] + x[i][j + 1] + x[i][j]);
        if (x_ <= spl) {
          uTmv.push_back(uL);
          vTmv.push_back(0.);
        } else {
          uTmv.push_back(uR);
          vTmv.push_back(0.);
        }
      }
      u.push_back(uTmv);
      v.push_back(vTmv);
    }
  };
  initializeConditions.emplace_back(initialCond);
}
void RiemannProblem1D::createEOSes(double gamma) {
  eoses.push_back(std::make_unique<EOSIdeal>(gamma));
}
