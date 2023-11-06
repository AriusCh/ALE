#include "problem.hpp"

#include <cassert>
#include <format>
#include <fstream>

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
      bottomBoundaryType(bottomBoundaryType_),
      writer(name_) {}

std::shared_ptr<GridALE> Problem::createALEGrid(int sizeX, int sizeY) const {
  return std::make_shared<GridALE>(sizeX, sizeY, xmin, xmax, ymin, ymax,
                                   uInitializer, vInitializer, rhoInitializer,
                                   pInitializer, eos);
}
void Problem::dumpGrid(std::shared_ptr<GridALE> grid, double t) const {
  std::string filename = std::format("{}_{:.3f}.txt", name, t);
  std::function<void(std::ofstream ofs)> outputF =
      [x = grid->x, y = grid->y, u = grid->u, v = grid->v, rho = grid->rho,
       p = grid->p, sizeX = grid->sizeX, sizeY = grid->sizeY,
       eos = grid->eos](std::ofstream ofs) {
        ofs << sizeX << " " << sizeY << std::endl;
        for (int i = 0; i < sizeX; i++) {
          for (int j = 0; j < sizeY; j++) {
            double xij =
                0.25 * (x[i][j] + x[i + 1][j] + x[i][j + 1] + x[i + 1][j + 1]);
            double yij =
                0.25 * (y[i][j] + y[i + 1][j] + y[i][j + 1] + y[i + 1][j + 1]);
            double uij =
                0.25 * (u[i][j] + u[i + 1][j] + u[i][j + 1] + u[i + 1][j + 1]);
            double vij =
                0.25 * (v[i][j] + v[i + 1][j] + v[i][j + 1] + v[i + 1][j + 1]);
            double rhoij = rho[i][j];
            double pij = p[i][j];
            double eij = eos->gete(rhoij, pij);
            auto output = std::format("{} {} {} {} {} {} {} {} {}", i, j, xij,
                                      yij, uij, vij, rhoij, pij, eij);
            ofs << output << std::endl;
          }
        }
      };
  writer.dumpData(std::move(outputF), filename);
}
void Problem::createProblem() {
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

  createProblem();
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
