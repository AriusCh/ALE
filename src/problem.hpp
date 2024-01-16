#ifndef ALE_SOLVER_SRC_PROBLEM_HPP_
#define ALE_SOLVER_SRC_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "boundary.hpp"
#include "eos.hpp"
#include "logger.hpp"
#include "output_mgr.hpp"

enum class ProblemDimension { e1D, e2D };

class Problem {
 public:
  Problem(
      const std::string &name_, double xmin_, double xmax_, double ymin_,
      double ymax_, double tmin_, double tmax_, BoundaryType leftBoundaryType_,
      BoundaryType topBoundaryType_, BoundaryType rightBoundaryType_,
      BoundaryType bottomBoundaryType_,
      std::function<double(double x, double y)> uInitializer_,
      std::function<double(double x, double y)> vInitializer_,
      std::function<double(double x, double y)> rhoInitializer_,
      std::function<double(double x, double y)> pInitializer_,
      std::function<std::shared_ptr<EOS>(double x, double y)> eosInitializer_,
      ProblemDimension dimension_);
  Problem(Problem const &rhs) = default;
  Problem(Problem &&rhs) = default;

  Problem &operator=(Problem const &rhs) = delete;
  Problem &operator=(Problem &&rhs) = delete;

  virtual ~Problem() = default;

 public:
 protected:
 private:
 public:
  const std::string name;  // Problem name

  const double xmin, xmax, ymin,
      ymax;                 // Starting rectangle simulation box boundaries
  const double tmin, tmax;  // Simulation time interval
  const BoundaryType leftBoundaryType;    // Left boundary type
  const BoundaryType topBoundaryType;     // Top boundary type
  const BoundaryType rightBoundaryType;   // Right boundary type
  const BoundaryType bottomBoundaryType;  // Bottom boundary type

  const std::function<double(double x, double y)>
      uInitializer;  // Function that returns u initial value depending on the
                     // coords
  const std::function<double(double x, double y)>
      vInitializer;  // Function that returns v initial value depending on the
                     // coords
  const std::function<double(double x, double y)>
      rhoInitializer;  // Function that returns rho initial values
  const std::function<double(double x, double y)>
      pInitializer;  // Function that returns p initial values
  const std::function<std::shared_ptr<EOS>(double x, double y)>
      eosInitializer;  // Function that returns EOSes to use in the grid

  const ProblemDimension dimension;

 protected:
  Logger logger;
};

/*  */
class RiemannProblem1Dx : public Problem {
 public:
  RiemannProblem1Dx(const std::string &name, double xmin, double xmax,
                    double tmax, double rhoL, double uL, double pL, double rhoR,
                    double uR, double pR, double spl, double gamma);

 public:
 protected:
 private:
};

// class CircularRiemannProblem : public Problem {
//  public:
//   CircularRiemannProblem(const std::string &name, double xmin, double xmax,
//                          double ymin, double ymax, double tmax, double rhoL,
//                          double uL, double vL, double pL, double rhoR,
//                          double uR, double vR, double pR, double spl,
//                          double gamma);
//
//  public:
//
//  private:
// };

#endif
