#include "problem.hpp"

#include <cassert>

Problem::Problem(ProblemType type_, AxisymmetryType symType_)
    : type(type_), symType(symType_) {}

std::unique_ptr<Grid> Problem::createGrid(size_t nx, size_t ny) const {
  std::unique_ptr<GridALE> grid = std::make_unique<GridALE>(nx, ny, calcRegion);

  for (int i = 0; i < internalAreas.size(); i++) {

  }
}

RiemannProblem1D::RiemannProblem1D(double xmin, double xmax, double tmax,
                                   double rhoL, double uL, double pL,
                                   double rhoR, double uR, double pR,
                                   double spl)
    : Problem(ProblemType::eRiemannProblem1D,
              AxisymmetryType::eNonSymmetrical) {
  assert(spl >= xmin && spl <= xmax);

  calcRegion =
      std::make_unique<PrimitiveRectangle>(xmin, 0.0, xmax - xmin, 1.0);
  calcRegion->setBoundaries(
      std::vector<BoundaryType>(4, BoundaryType::eExternalTransparent));

  internalAreas.push_back(
      std::make_unique<PrimitiveRectangle>(xmin, 0.0, spl - xmin, 1.0));

  internalAreas.push_back(
      std::make_unique<PrimitiveRectangle>(spl, 0.0, xmax - spl, 1.0));
}
