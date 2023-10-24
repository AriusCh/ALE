#ifndef ALE_SOLVER_SRC_PROBLEM_HPP_
#define ALE_SOLVER_SRC_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "eos.hpp"
#include "grid.hpp"
#include "logger.hpp"
#include "primitive.hpp"

enum class ProblemType { eRiemannProblem1D };

enum class AxisymmetryType { eNonSymmetrical, eSymmetrical };

/*  */
class Problem {
 public:
  Problem(ProblemType type_, AxisymmetryType symType_);
  Problem(Problem const &rhs) = delete;
  Problem(Problem &&rhs) = default;

  Problem &operator=(Problem const &rhs) = delete;
  Problem &operator=(Problem &&rhs) = default;

  virtual ~Problem() = 0;

 public:
  ProblemType getType() const;
  AxisymmetryType getSymmetryType() const;

  std::unique_ptr<Grid> createGrid(size_t nx, size_t ny) const;

 protected:
  virtual void createGeometry();
  virtual void createBoundaryTypes();
  virtual void createInitialConditions();
  virtual void createEOSes();

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
class RiemannProblem1D : Problem {
 public:
  RiemannProblem1D(double xmin, double xmax, double tmax, double rhoL,
                   double uL, double pL, double rhoR, double uR, double pR,
                   double spl, double gamma);

 protected:
  virtual void createGeometry(double xmin, double xmax);
  virtual void createBoundaryTypes() override;
  virtual void createInitialConditions(double rhoL, double pL, double uL,
                                       double rhoR, double pR, double uR,
                                       double spl);
  virtual void createEOSes(double gamma);
};

#endif
