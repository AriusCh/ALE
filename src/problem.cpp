#include "problem.hpp"

#include <cassert>
#include <format>
#include <fstream>

Problem::Problem(
    const std::string &name_, double xmin_, double xmax_, double ymin_,
    double ymax_, double tmin_, double tmax_, BoundaryType leftBoundaryType_,
    BoundaryType topBoundaryType_, BoundaryType rightBoundaryType_,
    BoundaryType bottomBoundaryType_,
    std::function<double(double x, double y)> uInitializer_,
    std::function<double(double x, double y)> vInitializer_,
    std::function<double(double x, double y)> rhoInitializer_,
    std::function<double(double x, double y)> pInitializer_,
    std::function<std::shared_ptr<EOS>(double x, double y)> eosInitializer_,
    ProblemDimension dimension_)
    : name(name_),
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
      uInitializer(uInitializer_),
      vInitializer(vInitializer_),
      rhoInitializer(rhoInitializer_),
      pInitializer(pInitializer_),
      eosInitializer(eosInitializer_),
      dimension(dimension_) {}

RiemannProblem1Dx::RiemannProblem1Dx(const std::string &name, double xmin,
                                     double xmax, double tmax, double rhoL,
                                     double uL, double pL, double rhoR,
                                     double uR, double pR, double spl,
                                     double gamma)
    : Problem(
          name, xmin, xmax, 0.0, 1.0, 0.0, tmax, BoundaryType::eWall,
          BoundaryType::eWall, BoundaryType::eWall, BoundaryType::eWall,
          [uL, uR, spl](double x, [[maybe_unused]] double y) {
            if (x <= spl) {
              return uL;
            } else {
              return uR;
            }
          },
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return 0.0;
          },
          [rhoL, rhoR, spl](double x, [[maybe_unused]] double y) {
            if (x <= spl) {
              return rhoL;

            } else {
              return rhoR;
            }
          },
          [pL, pR, spl](double x, [[maybe_unused]] double y) {
            if (x <= spl) {
              return pL;
            } else {
              return pR;
            }
          },
          [gamma]([[maybe_unused]] double x, [[maybe_unused]] double y) {
            static std::shared_ptr<EOSIdeal> eos =
                std::make_shared<EOSIdeal>(gamma);
            return eos;
          },
          ProblemDimension::e1D) {
  assert(spl >= xmin && spl <= xmax);
}
// void RiemannProblem1Dx::dumpGrid(std::shared_ptr<GridALE> grid,
//                                  double t) const {
//   std::string filename = std::format("{}_{:.3f}.txt", name, t);
//   std::function<void(std::ofstream ofs)> outputF =
//       [x = grid->x, y = grid->y, u = grid->u, v = grid->v, rho = grid->rho,
//        p = grid->p, sizeX = grid->sizeX, sizeY = grid->sizeY,
//        eos = grid->eos](std::ofstream ofs) {
//         ofs << sizeX << std::endl;
//         for (int i = 0; i < sizeX; i++) {
//           int j = sizeY / 2;
//           double xij =
//               0.25 * (x[i][j] + x[i + 1][j] + x[i][j + 1] + x[i + 1][j + 1]);
//           double uij =
//               0.25 * (u[i][j] + u[i + 1][j] + u[i][j + 1] + u[i + 1][j + 1]);
//           double vij =
//               0.25 * (v[i][j] + v[i + 1][j] + v[i][j + 1] + v[i + 1][j + 1]);
//           double rhoij = rho[i][j];
//           double pij = p[i][j];
//           double eij = eos->gete(rhoij, pij);
//           auto output =
//               std::format("{} {} {} {} {} {}", xij, uij, vij, rhoij, pij,
//               eij);
//           ofs << output << std::endl;
//         }
//       };
//   writer.dumpData(std::move(outputF), filename);
// }

CircularRiemannProblem::CircularRiemannProblem(
    const std::string &name, double xmin, double xmax, double ymin, double ymax,
    double tmax, double rhoL, double uL, double vL, double pL, double rhoR,
    double uR, double vR, double pR, double spl, double gamma)
    : Problem(
          name, xmin, xmax, ymin, ymax, 0.0, tmax, BoundaryType::eWall,
          BoundaryType::eWall, BoundaryType::eWall, BoundaryType::eWall,
          [uL, uR, spl](double x, double y) {
            if (x * x + y * y <= spl * spl) {
              return uL;
            } else {
              return uR;
            }
          },
          [vL, vR, spl](double x, double y) {
            if (x * x + y * y <= spl * spl) {
              return vL;
            } else {
              return vR;
            }
          },
          [rhoL, rhoR, spl](double x, double y) {
            if (x * x + y * y <= spl * spl) {
              return rhoL;
            } else {
              return rhoR;
            }
          },
          [pL, pR, spl](double x, double y) {
            if (x * x + y * y <= spl * spl) {
              return pL;
            } else {
              return pR;
            }
          },
          [gamma]([[maybe_unused]] double x, [[maybe_unused]] double y) {
            static std::shared_ptr<EOSIdeal> eos =
                std::make_shared<EOSIdeal>(gamma);
            return eos;
          },
          ProblemDimension::e1D) {
  assert(spl >= xmin && spl <= xmax);
}

// Riemann test problems
std::shared_ptr<Problem> sodTest = std::make_shared<RiemannProblem1Dx>(
    "sod-test", 0.0, 1.0, 0.2, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 0.5, 1.4);
std::shared_ptr<Problem> blastWave2D = std::make_shared<CircularRiemannProblem>(
    "blastWave", -1.0, 1.0, -1.0, 1.0, 0.25, 1.0, 0.0, 0.0, 1.0, 0.125, 0.0,
    0.0, 0.1, 0.4, 1.4);
std::shared_ptr<Problem> toro1x = std::make_shared<RiemannProblem1Dx>(
    "toro1x", 0.0, 1.0, 0.2, 1.0, 0.75, 1.0, 0.125, 0.0, 0.1, 0.3, 1.4);
std::shared_ptr<Problem> toro2x = std::make_shared<RiemannProblem1Dx>(
    "toro2x", 0.0, 1.0, 0.15, 1.0, -2.0, 0.4, 1.0, 2.0, 0.4, 0.5, 1.4);
std::shared_ptr<Problem> toro3x = std::make_shared<RiemannProblem1Dx>(
    "toro3x", 0.0, 1.0, 0.012, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.01, 0.5, 1.4);
std::shared_ptr<Problem> toro4x = std::make_shared<RiemannProblem1Dx>(
    "toro4x", 0.0, 1.0, 0.035, 5.99924, 19.5975, 460.894, 5.99242, -6.19633,
    46.0950, 0.4, 1.4);
std::shared_ptr<Problem> toro5x = std::make_shared<RiemannProblem1Dx>(
    "toro5x", 0.0, 1.0, 0.012, 1.0, -19.59745, 1000.0, 1.0, -19.59745, 0.01,
    0.8, 1.4);
