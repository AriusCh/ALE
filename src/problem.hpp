#ifndef ALE_SOLVER_SRC_PROBLEM_HPP_
#define ALE_SOLVER_SRC_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "boundary.hpp"
#include "eos.hpp"
#include "grid.hpp"
#include "grid_manager.hpp"
#include "logger.hpp"
#include "primitive.hpp"

enum class ProblemType { eRiemannProblem1Dx };

enum class AxisymmetryType { eNonSymmetrical, eSymmetrical };

/*  */
class Problem {
 public:
  Problem(ProblemType type_, AxisymmetryType symType_, double tmin,
          double tmax);
  Problem(Problem const &rhs) = delete;
  Problem(Problem &&rhs) = default;

  Problem &operator=(Problem const &rhs) = delete;
  Problem &operator=(Problem &&rhs) = default;

  virtual ~Problem() = 0;

 public:
  ProblemType getType() const;
  AxisymmetryType getSymmetryType() const;

  std::vector<std::unique_ptr<GridALE>> createALEGrids(size_t nx, size_t ny);

  std::vector<BoundaryType> getLeftBoundaryTypes() const;
  std::vector<BoundaryType> getTopBoundaryTypes() const;
  std::vector<BoundaryType> getRightBoundaryTypes() const;
  std::vector<BoundaryType> getBottomBoundaryTypes() const;

  double getTMin() const;
  double getTMax() const;

 protected:
  void create();
  virtual void createGeometry() = 0;       // create region of calculation
  virtual void createBoundaryTypes() = 0;  // create boundary types for each
                                           // polygon in region
  virtual void
  createInitialConditions() = 0;  // create initial condition functions for each
                                  // polygon in region
  virtual void
  createEOSes() = 0;  // create equation of state for each polygon in region

 private:
 protected:
  std::unique_ptr<Region2D> calcRegion;
  std::vector<BoundaryType> leftBoundaryTypes;
  std::vector<BoundaryType> topBoundaryTypes;
  std::vector<BoundaryType> rightBoundaryTypes;
  std::vector<BoundaryType> bottomBoundaryTypes;
  std::vector<std::function<void(
      const std::vector<std::vector<double>> &,
      const std::vector<std::vector<double>> &,
      std::vector<std::vector<double>> &, std::vector<std::vector<double>> &,
      std::vector<std::vector<double>> &, std::vector<std::vector<double>> &)>>
      initializeConditions;
  std::vector<std::unique_ptr<EOS>> eoses;

  double tmin = 0., tmax;

  Logger logger;

 private:
  ProblemType type;
  AxisymmetryType symType;
};

/*  */
class RiemannProblem1Dx : public Problem {
 public:
  RiemannProblem1Dx(double xmin, double xmax, double tmax, double rhoL,
                    double uL, double pL, double rhoR, double uR, double pR,
                    double spl, double gamma);

 protected:
  virtual void createGeometry() override;
  virtual void createBoundaryTypes() override;
  virtual void createInitialConditions() override;
  virtual void createEOSes() override;

 private:
  double xmin, xmax;
  double rhoL, pL, uL;
  double rhoR, pR, uR;
  double spl;
  double gamma;
};

#endif
