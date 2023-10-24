#ifndef ALE_SOLVER_SRC_PROBLEM_HPP_
#define ALE_SOLVER_SRC_PROBLEM_HPP_

#include <functional>
#include <memory>
#include <vector>

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
  std::unique_ptr<Region2D> calcRegion;
  std::vector<std::function<void(
      const std::vector<std::vector<double>> &,
      const std::vector<std::vector<double>> &,
      std::vector<std::vector<double>> &, std::vector<std::vector<double>> &,
      std::vector<std::vector<double>> &, std::vector<std::vector<double>> &)>>
      initialConditions;

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
                   double spl);
};

#endif
