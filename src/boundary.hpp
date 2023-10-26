#ifndef ALE_SOLVER_SRC_BOUNDARY_HPP_
#define ALE_SOLVER_SRC_BOUNDARY_HPP_

#include <vector>

enum class BoundaryType {
  eExternalTransparent,
  eExternalNoSlipWall,
  eExternalFree
};

enum class ExternalBoundarySide { eLeft, eTop, eRight, eBottom };

class Boundary {
 public:
  Boundary(BoundaryType type_);

 public:
  BoundaryType getType() const;

 private:
  BoundaryType type;
};

class ExternalBoundary : public Boundary {
 public:
  ExternalBoundary(BoundaryType type_, ExternalBoundarySide side_);

 protected:
  ExternalBoundarySide side;
  std::vector<std::vector<double>> &x;
  std::vector<std::vector<double>> &y;
  std::vector<std::vector<double>> &rho;
  std::vector<std::vector<double>> &p;
  std::vector<std::vector<double>> &m;
};

#endif
