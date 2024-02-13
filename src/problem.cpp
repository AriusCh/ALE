#include "problem.hpp"

#include <cassert>
#include <format>
#include <fstream>

Problem::Problem(
    const std::string &name_, double xmin_, double xmax_, double ymin_,
    double ymax_, double tmin_, double tmax_, const std::deque<double> &tOut_,
    double tMul_, BoundaryType leftBoundaryType_, BoundaryType topBoundaryType_,
    BoundaryType rightBoundaryType_, BoundaryType bottomBoundaryType_,
    ProblemDimension dimension_,
    std::function<double(double x, double y)> uInitializer_,
    std::function<double(double x, double y)> vInitializer_,
    std::function<double(double x, double y)> rhoInitializer_,
    std::function<double(double x, double y)> pInitializer_,
    std::function<std::shared_ptr<EOS>(double x, double y)> eosInitializer_)
    : name(name_),
      xmin(xmin_),
      xmax(xmax_),
      ymin(ymin_),
      ymax(ymax_),
      tmin(tmin_),
      tmax(tmax_),
      tOut(tOut_),
      tMul(tMul_),
      leftBoundaryType(leftBoundaryType_),
      topBoundaryType(topBoundaryType_),
      rightBoundaryType(rightBoundaryType_),
      bottomBoundaryType(bottomBoundaryType_),
      dimension(dimension_),
      uInitializer(uInitializer_),
      vInitializer(vInitializer_),
      rhoInitializer(rhoInitializer_),
      pInitializer(pInitializer_),
      eosInitializer(eosInitializer_) {}

RiemannProblem1Dx::RiemannProblem1Dx(const std::string &name, double xmin,
                                     double xmax, double tmax,
                                     const std::deque<double> &tOut,
                                     double rhoL, double uL, double pL,
                                     double rhoR, double uR, double pR,
                                     double spl, double gamma)
    : Problem(
          name, xmin, xmax, 0.0, 1.0, 0.0, tmax, tOut, 1.0, BoundaryType::eWall,
          BoundaryType::eWall, BoundaryType::eWall, BoundaryType::eWall,
          ProblemDimension::e1D,
          [uL, uR, spl](double x, [[maybe_unused]] double y) {
            if (x < spl) {
              return uL;
            } else {
              return uR;
            }
          },
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return 0.0;
          },
          [rhoL, rhoR, spl](double x, [[maybe_unused]] double y) {
            if (x < spl) {
              return rhoL;

            } else {
              return rhoR;
            }
          },
          [pL, pR, spl](double x, [[maybe_unused]] double y) {
            if (x < spl) {
              return pL;
            } else {
              return pR;
            }
          },
          [gamma]([[maybe_unused]] double x, [[maybe_unused]] double y) {
            static std::shared_ptr<EOSIdeal> eos =
                std::make_shared<EOSIdeal>(gamma);
            return eos;
          }) {
  assert(spl >= xmin && spl <= xmax);
}

CircularRiemannProblem::CircularRiemannProblem(
    const std::string &name, double xmin, double xmax, double ymin, double ymax,
    double tmax, const std::deque<double> &tOut, double rhoL, double uL,
    double vL, double pL, double rhoR, double uR, double vR, double pR,
    double spl, double gamma)
    : Problem(
          name, xmin, xmax, ymin, ymax, 0.0, tmax, tOut, 1.0,
          BoundaryType::eWall, BoundaryType::eWall, BoundaryType::eWall,
          BoundaryType::eWall, ProblemDimension::e2D,
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
          }) {
  assert(spl >= xmin && spl <= xmax);
}

LaserVolumeTargetProblem::LaserVolumeTargetProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::deque<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin)
    : Problem(
          name, xmin, xmax, ymin, 0.0, 0.0, tmax, tOut, 1.0e12,
          BoundaryType::eWall, BoundaryType::eFree, BoundaryType::eWall,
          BoundaryType::eWall, ProblemDimension::e2D,
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return 0.0;
          },
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return 0.0;
          },
          [rhoM]([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return rhoM;
          },
          [pCold, pHeat, RL, dSkin](double x, double y) {
            if (std::abs(x) <= RL && y > -dSkin) {
              return pHeat;
            } else {
              return pCold;
            }
          },
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            static std::shared_ptr<EOSMGAlPrecise6> eos =
                std::make_shared<EOSMGAlPrecise6>();
            return eos;
          }) {}

TriplePointShock::TriplePointShock(const std::string &name, double xmin,
                                   double xmax, double ymin, double ymax,
                                   double tmax, const std::deque<double> &tOut,
                                   double xLeft, double yTop, double rhoLeft,
                                   double pLeft, double rhoBottom,
                                   double pBottom, double rhoTop, double pTop,
                                   double gammaLeft, double gammaBottom,
                                   double gammaTop)
    : Problem(
          name, xmin, xmax, ymin, ymax, 0.0, tmax, tOut, 1.0,
          BoundaryType::eWall, BoundaryType::eWall, BoundaryType::eWall,
          BoundaryType::eWall, ProblemDimension::e2D,
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return 0.0;
          },
          []([[maybe_unused]] double x, [[maybe_unused]] double y) {
            return 0.0;
          },
          [rhoLeft, rhoTop, rhoBottom, xLeft, yTop](double x, double y) {
            if (x <= xLeft) {
              return rhoLeft;
            }
            if (y <= yTop) {
              return rhoBottom;
            }
            return rhoTop;
          },
          [pLeft, pTop, pBottom, xLeft, yTop](double x, double y) {
            if (x <= xLeft) {
              return pLeft;
            }
            if (y <= yTop) {
              return pBottom;
            }
            return pTop;
          },

          [gammaLeft, gammaTop, gammaBottom, xLeft, yTop](double x, double y) {
            static std::shared_ptr<EOSIdeal> eosLeft =
                std::make_shared<EOSIdeal>(gammaLeft);
            static std::shared_ptr<EOSIdeal> eosBottom =
                std::make_shared<EOSIdeal>(gammaBottom);
            static std::shared_ptr<EOSIdeal> eosTop =
                std::make_shared<EOSIdeal>(gammaTop);
            if (x <= xLeft) {
              return eosLeft;
            }
            if (y <= yTop) {
              return eosBottom;
            }
            return eosTop;
          }

      ) {
  assert(xLeft > xmin && xLeft < xmax);
  assert(yTop > ymin && yTop < ymax);
}

// Riemann test problems
std::shared_ptr<Problem> sodTest = std::make_shared<RiemannProblem1Dx>(
    "sod-test", 0.0, 1.0, 0.2, std::deque<double>{}, 1.0, 0.0, 1.0, 0.125, 0.0,
    0.1, 0.5, 5.0 / 3.0);
std::shared_ptr<Problem> blastWave2D = std::make_shared<CircularRiemannProblem>(
    "blastWave", 0.0, 1.0, 0.0, 1.0, 0.25, std::deque<double>{}, 1.0, 0.0, 0.0,
    1.0, 0.125, 0.0, 0.0, 0.1, 0.4, 1.4);
std::shared_ptr<Problem> toro1x = std::make_shared<RiemannProblem1Dx>(
    "toro1x", 0.0, 1.0, 0.2, std::deque<double>{}, 1.0, 0.75, 1.0, 0.125, 0.0,
    0.1, 0.3, 1.4);
std::shared_ptr<Problem> toro2x = std::make_shared<RiemannProblem1Dx>(
    "toro2x", 0.0, 1.0, 0.15, std::deque<double>{}, 1.0, -2.0, 0.4, 1.0, 2.0,
    0.4, 0.5, 1.4);
std::shared_ptr<Problem> toro3x = std::make_shared<RiemannProblem1Dx>(
    "toro3x", 0.0, 1.0, 0.012, std::deque<double>{}, 1.0, 0.0, 1000.0, 1.0, 0.0,
    0.01, 0.5, 1.4);
std::shared_ptr<Problem> toro4x = std::make_shared<RiemannProblem1Dx>(
    "toro4x", 0.0, 1.0, 0.035, std::deque<double>{}, 5.99924, 19.5975, 460.894,
    5.99242, -6.19633, 46.0950, 0.4, 1.4);
std::shared_ptr<Problem> toro5x = std::make_shared<RiemannProblem1Dx>(
    "toro5x", 0.0, 1.0, 0.012, std::deque<double>{}, 1.0, -19.59745, 1000.0,
    1.0, -19.59745, 0.01, 0.8, 1.4);

// LaserVolumeDefaultProblem
std::shared_ptr<Problem> laserVolumeTargetProblem =
    std::make_shared<LaserVolumeTargetProblem>(
        "laser-al", 0.0, 1280.e-9, -800.e-9, 115.2e-12,
        std::deque<double>{0.1 * 9.6e-12, 1.0 * 9.6e-12, 2.0 * 9.6e-12,
                           3.0 * 9.6e-12, 4.0 * 9.6e-12, 5.0 * 9.6e-12,
                           6.0 * 9.6e-12, 7.0 * 9.6e-12, 8.0 * 9.6e-12,
                           9.0 * 9.6e-12, 10.0 * 9.6e-12, 11.0 * 9.6e-12},
        2413., 0, 35.6e9, 200e-9, 80e-9);

// Triple point shock wave
std::shared_ptr<Problem> triplePointShock = std::make_shared<TriplePointShock>(
    "triple-point", 0.0, 7.0, 0.0, 3.0, 5.0,
    std::deque<double>{0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.3, 3.5, 4.0, 4.5}, 1.0,
    1.5, 1.0, 1.0, 1.0, 0.1, 0.125, 0.1, 1.5, 1.4, 1.5);
