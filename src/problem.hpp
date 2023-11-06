#ifndef ALE_SOLVER_SRC_PROBLEM_HPP_
#define ALE_SOLVER_SRC_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "boundary.hpp"
#include "eos.hpp"
#include "grid.hpp"
#include "logger.hpp"
#include "output_mgr.hpp"

enum class ProblemType { eRiemannProblem1Dx };

enum class AxisymmetryType { eNone, eCylindrical };

/*  */
class Problem {
 public:
  Problem(const std::string &name_, ProblemType type_, AxisymmetryType symType_,
          double xmin_, double xmax_, double ymin_, double ymax_, double tmin_,
          double tmax_, BoundaryType leftBoundaryType_,
          BoundaryType topBoundaryType_, BoundaryType rightBoundaryType_,
          BoundaryType bottomBoundaryType_);
  Problem(Problem const &rhs) = default;
  Problem(Problem &&rhs) = default;

  Problem &operator=(Problem const &rhs) = delete;
  Problem &operator=(Problem &&rhs) = delete;

  virtual ~Problem() = default;

 public:
  virtual std::shared_ptr<GridALE> createALEGrid(
      int sizeX,
      int sizeY) const;  // Create a rectangle ALE grid
  virtual void dumpGrid(std::shared_ptr<GridALE> grid,
                        double t) const;  // Dump grid as 2D

 protected:
  void createProblem();  // Function that calls all create functions
  virtual void createInitializers() = 0;  // create initial condition function
  virtual void createEOS() = 0;           // create equation of state

 private:
 public:
  const std::string name;         // Problem name
  const ProblemType type;         // Problem type
  const AxisymmetryType symType;  // Problem axisymmetry

  const double xmin, xmax, ymin,
      ymax;                 // Starting rectangle simulation box boundaries
  const double tmin, tmax;  // Simulation time interval
  const BoundaryType leftBoundaryType;    // Left boundary type
  const BoundaryType topBoundaryType;     // Top boundary type
  const BoundaryType rightBoundaryType;   // Right boundary type
  const BoundaryType bottomBoundaryType;  // Bottom boundary type

  std::function<double(double x, double y)>
      uInitializer;  // Function that returns u initial value depending on the
                     // coords
  std::function<double(double x, double y)>
      vInitializer;  // Function that returns v initial value depending on the
                     // coords
  std::function<double(double x, double y)>
      rhoInitializer;  // Function that returns rho initial value
  std::function<double(double x, double y)>
      pInitializer;          // Function that returns p initial value
  std::shared_ptr<EOS> eos;  // EOS to use in the grid

 protected:
  Writer writer;
  Logger logger;
};

/*  */
class RiemannProblem1Dx : public Problem {
 public:
  RiemannProblem1Dx(const std::string &name, double xmin, double xmax,
                    double tmax, double rhoL, double uL, double pL, double rhoR,
                    double uR, double pR, double spl, double gamma);

 protected:
  virtual void createInitializers() override;
  virtual void createEOS() override;

 private:
  double rhoL, pL, uL;
  double rhoR, pR, uR;
  double spl;
  double gamma;
};

#endif
