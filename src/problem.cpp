#include "problem.hpp"

#include <cassert>

Problem::Problem(ProblemType type_, AxisymmetryType symType_, double tmin,
                 double tmax)
    : type(type_), symType(symType_), tmin(tmin), tmax(tmax) {}

ProblemType Problem::getType() const { return type; }
AxisymmetryType Problem::getSymmetryType() const { return symType; }

std::vector<std::unique_ptr<GridALE>> Problem::createALEGrids(size_t nx,
                                                              size_t ny) {
  const auto &polygons = calcRegion->getPolygons();

  std::vector<std::unique_ptr<GridALE>> grids(polygons.size());
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
    grids[i] =
        std::make_unique<GridALE>(x, y, rho, p, u, v, std::move(eoses[i]));
  }
  return grids;
}
void Problem::create() {
  createGeometry();
  createBoundaryTypes();
  createInitialConditions();
  createEOSes();
}

RiemannProblem1Dx::RiemannProblem1Dx(double xmin, double xmax, double tmax,
                                     double rhoL, double uL, double pL,
                                     double rhoR, double uR, double pR,
                                     double spl, double gamma)
    : xmin(xmin),
      xmax(xmax),
      rhoL(rhoL),
      pL(pL),
      uL(uL),
      rhoR(rhoR),
      pR(pR),
      uR(uR),
      Problem(ProblemType::eRiemannProblem1Dx, AxisymmetryType::eNonSymmetrical,
              0.0, tmax) {
  assert(spl >= xmin && spl <= xmax);

  create();
}
void RiemannProblem1Dx::createGeometry() {
  std::unique_ptr<Polygon> rect =
      std::make_unique<Rectangle>(xmin, 0.0, xmax - xmin, 1.0);
  calcRegion = std::make_unique<Region2D>(std::move(rect));
}
void RiemannProblem1Dx::createBoundaryTypes() {
  leftBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
  topBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
  rightBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
  bottomBoundaryTypes.push_back(BoundaryType::eExternalTransparent);
}
void RiemannProblem1Dx::createInitialConditions() {
  auto initialCond = [rhoL = this->rhoL, pL = this->pL, uL = this->uL,
                      rhoR = this->rhoR, pR = this->pR, uR = this->uR,
                      spl = this->spl](
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
void RiemannProblem1Dx::createEOSes() {
  eoses.push_back(std::make_unique<EOSIdeal>(gamma));
}
