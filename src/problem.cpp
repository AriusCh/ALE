#include "problem.hpp"

#include <cassert>
#include <format>
#include <fstream>

Problem Problem::createRiemannProblem1Dx(const std::string &name, double xmin,
                                         double xmax, double tmax,
                                         const std::deque<double> &tOut,
                                         double uL, double rhoL, double pL,
                                         double uR, double rhoR, double pR,
                                         double spl, double gamma) {
  assert(spl >= xmin && spl <= xmax);
  double ymin = 0.0;
  double ymax = 1.0;
  double tmin = 0.0;
  double tMul = 1.0;
  BoundaryType leftBoundaryType = BoundaryType::eWall;
  BoundaryType topBoundaryType = BoundaryType::eWall;
  BoundaryType rightBoundaryType = BoundaryType::eWall;
  BoundaryType bottomBoundaryType = BoundaryType::eWall;
  ProblemDimension dimension = ProblemDimension::e1D;
  auto uInitializer = [uL, uR, spl](double x, [[maybe_unused]] double y) {
    if (x < spl) {
      return uL;
    } else {
      return uR;
    }
  };
  auto vInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  auto rhoInitializer = [rhoL, rhoR, spl](double x, [[maybe_unused]] double y) {
    if (x < spl) {
      return rhoL;

    } else {
      return rhoR;
    }
  };
  auto pInitializer = [pL, pR, spl](double x, [[maybe_unused]] double y) {
    if (x < spl) {
      return pL;
    } else {
      return pR;
    }
  };
  auto eosInitializer = [gamma]([[maybe_unused]] double x,
                                [[maybe_unused]] double y) {
    static std::shared_ptr<EOSIdeal> eos = std::make_shared<EOSIdeal>(gamma);
    return eos;
  };
  return Problem(name, xmin, xmax, ymin, ymax, tmin, tmax, tOut, tMul,
                 leftBoundaryType, topBoundaryType, rightBoundaryType,
                 bottomBoundaryType, dimension, uInitializer, vInitializer,
                 rhoInitializer, pInitializer, eosInitializer);
}

Problem Problem::createCircularRiemannProblem(
    const std::string &name, double xmin, double xmax, double ymin, double ymax,
    double tmax, const std::deque<double> &tOut, double uIn, double vIn,
    double rhoIn, double pIn, double uOut, double vOut, double rhoOut,
    double pOut, double spl, double gamma) {
  assert(spl >= xmin && spl <= xmax);
  assert(spl >= ymin && spl <= ymax);

  double tmin = 0.0;
  double tMul = 1.0;
  BoundaryType leftBoundaryType = BoundaryType::eWall;
  BoundaryType topBoundaryType = BoundaryType::eWall;
  BoundaryType rightBoundaryType = BoundaryType::eWall;
  BoundaryType bottomBoundaryType = BoundaryType::eWall;
  ProblemDimension dimension = ProblemDimension::e2D;
  auto uInitializer = [uIn, uOut, spl](double x, double y) {
    if (x * x + y * y <= spl * spl) {
      return uIn;
    } else {
      return uOut;
    }
  };
  auto vInitializer = [vIn, vOut, spl](double x, double y) {
    if (x * x + y * y <= spl * spl) {
      return vIn;
    } else {
      return vOut;
    }
  };
  auto rhoInitializer = [rhoIn, rhoOut, spl](double x, double y) {
    if (x * x + y * y <= spl * spl) {
      return rhoIn;
    } else {
      return rhoOut;
    }
  };
  auto pInitializer = [pIn, pOut, spl](double x, double y) {
    if (x * x + y * y <= spl * spl) {
      return pIn;
    } else {
      return pOut;
    }
  };
  auto eosInitializer = [gamma]([[maybe_unused]] double x,
                                [[maybe_unused]] double y) {
    static std::shared_ptr<EOSIdeal> eos = std::make_shared<EOSIdeal>(gamma);
    return eos;
  };

  return Problem(name, xmin, xmax, ymin, ymax, tmin, tmax, tOut, tMul,
                 leftBoundaryType, topBoundaryType, rightBoundaryType,
                 bottomBoundaryType, dimension, uInitializer, vInitializer,
                 rhoInitializer, pInitializer, eosInitializer);
}
Problem Problem::createLaserVolumeTargetProblem(
    const std::string &name, double xmin, double xmax, double ymin, double tmax,
    const std::deque<double> &tOut, double rhoM, double pCold, double pHeat,
    double RL, double dSkin) {
  assert(ymin < -dSkin);

  double ymax = 0.0;
  double tmin = 0.0;
  double tMul = 1e12;
  BoundaryType leftBoundaryType = BoundaryType::eWall;
  BoundaryType topBoundaryType = BoundaryType::eFree;
  BoundaryType rightBoundaryType = BoundaryType::eWall;
  BoundaryType bottomBoundaryType = BoundaryType::eWall;
  ProblemDimension dimension = ProblemDimension::e2D;

  auto uInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  auto vInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  auto rhoInitializer = [rhoM]([[maybe_unused]] double x,
                               [[maybe_unused]] double y) { return rhoM; };
  auto pInitializer = [pCold, pHeat, RL, dSkin](double x, double y) {
    if (std::abs(x) <= RL && y > -dSkin) {
      return pHeat;
    } else {
      return pCold;
    }
  };
  auto eosInitializer = []([[maybe_unused]] double x,
                           [[maybe_unused]] double y) {
    static std::shared_ptr<EOSMGAlPrecise6> eos =
        std::make_shared<EOSMGAlPrecise6>();
    return eos;
  };

  return Problem(name, xmin, xmax, ymin, ymax, tmin, tmax, tOut, tMul,
                 leftBoundaryType, topBoundaryType, rightBoundaryType,
                 bottomBoundaryType, dimension, uInitializer, vInitializer,
                 rhoInitializer, pInitializer, eosInitializer);
}

Problem Problem::createTriplePointProblem(
    const std::string &name, double xmin, double xmax, double ymin, double ymax,
    double tmax, const std::deque<double> &tOut, double vSplit, double hSplit,
    double rhoLeft, double pLeft, double gammaLeft, double rhoTop, double pTop,
    double gammaTop, double rhoBottom, double pBottom, double gammaBottom) {
  assert(vSplit > xmin && vSplit < xmax);
  assert(hSplit > ymin && hSplit < ymax);

  double tmin = 0.0;
  double tMul = 1.0;

  BoundaryType leftBoundaryType = BoundaryType::eWall;
  BoundaryType topBoundaryType = BoundaryType::eWall;
  BoundaryType rightBoundaryType = BoundaryType::eWall;
  BoundaryType bottomBoundaryType = BoundaryType::eWall;
  ProblemDimension dimension = ProblemDimension::e2D;

  auto uInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  auto vInitializer = []([[maybe_unused]] double x, [[maybe_unused]] double y) {
    return 0.0;
  };
  auto rhoInitializer = [rhoLeft, rhoTop, rhoBottom, vSplit, hSplit](double x,
                                                                     double y) {
    if (x <= vSplit) {
      return rhoLeft;
    }
    if (y <= hSplit) {
      return rhoBottom;
    }
    return rhoTop;
  };
  auto pInitializer = [pLeft, pTop, pBottom, vSplit, hSplit](double x,
                                                             double y) {
    if (x <= vSplit) {
      return pLeft;
    }
    if (y <= hSplit) {
      return pBottom;
    }
    return pTop;
  };

  auto eosInitializer = [gammaLeft, gammaTop, gammaBottom, vSplit, hSplit](
                            double x, double y) {
    static std::shared_ptr<EOSIdeal> eosLeft =
        std::make_shared<EOSIdeal>(gammaLeft);
    static std::shared_ptr<EOSIdeal> eosBottom =
        std::make_shared<EOSIdeal>(gammaBottom);
    static std::shared_ptr<EOSIdeal> eosTop =
        std::make_shared<EOSIdeal>(gammaTop);
    if (x <= vSplit) {
      return eosLeft;
    }
    if (y <= hSplit) {
      return eosBottom;
    }
    return eosTop;
  };

  return Problem(name, xmin, xmax, ymin, ymax, tmin, tmax, tOut, tMul,
                 leftBoundaryType, topBoundaryType, rightBoundaryType,
                 bottomBoundaryType, dimension, uInitializer, vInitializer,
                 rhoInitializer, pInitializer, eosInitializer);
}

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

// std::shared_ptr<Problem> toro1x = std::make_shared<RiemannProblem1Dx>(
//     "toro1x", 0.0, 1.0, 0.2, std::deque<double>{}, 1.0, 0.75, 1.0, 0.125,
//     0.0, 0.1, 0.3, 1.4);
// std::shared_ptr<Problem> toro2x = std::make_shared<RiemannProblem1Dx>(
//     "toro2x", 0.0, 1.0, 0.15, std::deque<double>{}, 1.0, -2.0, 0.4, 1.0, 2.0,
//     0.4, 0.5, 1.4);
// std::shared_ptr<Problem> toro3x = std::make_shared<RiemannProblem1Dx>(
//     "toro3x", 0.0, 1.0, 0.012, std::deque<double>{}, 1.0, 0.0, 1000.0, 1.0,
//     0.0, 0.01, 0.5, 1.4);
// std::shared_ptr<Problem> toro4x = std::make_shared<RiemannProblem1Dx>(
//     "toro4x", 0.0, 1.0, 0.035, std::deque<double>{}, 5.99924, 19.5975,
//     460.894, 5.99242, -6.19633, 46.0950, 0.4, 1.4);
// std::shared_ptr<Problem> toro5x = std::make_shared<RiemannProblem1Dx>(
//     "toro5x", 0.0, 1.0, 0.012, std::deque<double>{}, 1.0, -19.59745, 1000.0,
//     1.0, -19.59745, 0.01, 0.8, 1.4);

// Riemann test problems
Problem Problems::sodTest = Problem::createRiemannProblem1Dx(
    "sod-test", 0.0, 1.0, 0.2, std::deque<double>{}, 0.0, 1.0, 1.0, 0.0, 0.125,
    0.1, 0.5, 5.0 / 3.0);

Problem Problems::blastWave2D = Problem::createCircularRiemannProblem(
    "blastWave", 0.0, 1.0, 0.0, 1.0, 0.25, std::deque<double>{}, 0.0, 0.0, 1.0,
    1.0, 0.0, 0.0, 0.125, 0.1, 0.4, 1.4);

// LaserVolumeDefaultProblem
Problem Problems::laserVolumeTarget = Problem::createLaserVolumeTargetProblem(
    "laser-al", 0.0, 1280e-9, -800e-9, 115.2e-12,
    std::deque<double>{0.1 * 9.6e-12, 1.0 * 9.6e-12, 2.0 * 9.6e-12,
                       3.0 * 9.6e-12, 4.0 * 9.6e-12, 5.0 * 9.6e-12,
                       6.0 * 9.6e-12, 7.0 * 9.6e-12, 8.0 * 9.6e-12,
                       9.0 * 9.6e-12, 10.0 * 9.6e-12, 11.0 * 9.6e-12},
    2413.0, 0.0, 35.6e9, 200e-9, 80e-9);

// Triple point shock wave
Problem Problems::triplePointShock = Problem::createTriplePointProblem(
    "triple-point", 0.0, 7.0, 0.0, 3.0, 5.0,
    std::deque<double>{0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.3, 3.5, 4.0, 4.5}, 1.0,
    1.5, 1.0, 1.0, 1.5, 0.125, 0.1, 1.5, 1.0, 0.1, 1.4);
