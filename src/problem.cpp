#include "problem.hpp"

#include <cassert>

Problem::Problem(const std::string &name_, ProblemType type_,
                 AxisymmetryType symType_, double xmin_, double xmax_,
                 double ymin_, double ymax_, double tmin_, double tmax_,
                 BoundaryType leftBoundaryType_, BoundaryType topBoundaryType_,
                 BoundaryType rightBoundaryType_,
                 BoundaryType bottomBoundaryType_)
    : name(name_),
      type(type_),
      symType(symType_),
      xmin(xmin_),
      xmax(xmax_),
      ymin(ymin_),
      ymax(ymax_),
      tmin(tmin_),
      tmax(tmax_),
      leftBoundaryType(leftBoundaryType_),
      topBoundaryType(topBoundaryType_),
      rightBoundaryType(rightBoundaryType_),
      bottomBoundaryType(bottomBoundaryType_) {}

std::shared_ptr<GridALE> Problem::createALEGrid(int sizeX, int sizeY) const {
  return std::make_shared<GridALE>(sizeX, sizeY, xmin, xmax, ymin, ymax,
                                   uInitializer, vInitializer, rhoInitializer,
                                   pInitializer, eos);
}
void Problem::create() {
  createInitializers();
  createEOS();
}

RiemannProblem1Dx::RiemannProblem1Dx(const std::string &name, double xmin,
                                     double xmax, double tmax, double rhoL,
                                     double uL, double pL, double rhoR,
                                     double uR, double pR, double spl,
                                     double gamma)
    : Problem(name, ProblemType::eRiemannProblem1Dx, AxisymmetryType::eNone,
              xmin, xmax, 0.0, 1.0, 0.0, tmax,
              BoundaryType::eExternalTransparent,
              BoundaryType::eExternalTransparent,
              BoundaryType::eExternalTransparent,
              BoundaryType::eExternalTransparent),
      rhoL(rhoL),
      pL(pL),
      uL(uL),
      rhoR(rhoR),
      pR(pR),
      uR(uR),
      spl(spl),
      gamma(gamma) {
  assert(spl >= xmin && spl <= xmax);

  create();
}
void RiemannProblem1Dx::createInitializers() {
  uInitializer = [uL = this->uL, uR = this->uR, spl = this->spl](
                     double x, [[maybe_unused]] double y) {
    if (x <= spl) {
      return uL;
    } else {
      return uR;
    }
  };
  vInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  rhoInitializer = [rhoL = this->rhoL, rhoR = this->rhoR, spl = this->spl](
                       double x, [[maybe_unused]] double y) {
    if (x <= spl) {
      return rhoL;
    } else {
      return rhoR;
    }
  };
  pInitializer = [pL = this->pL, pR = this->pR, spl = this->spl](
                     double x, [[maybe_unused]] double y) {
    if (x <= spl) {
      return pL;
    } else {
      return pR;
    }
  };
}
void RiemannProblem1Dx::createEOS() { eos = std::make_shared<EOSIdeal>(gamma); }

// Riemann test problems
std::shared_ptr<Problem> sodTest = std::make_shared<RiemannProblem1Dx>(
    "sod-test", 0.0, 1.0, 0.2, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 0.5, 1.4);
