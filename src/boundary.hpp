#ifndef ALE_SOLVER_SRC_BOUNDARY_HPP_
#define ALE_SOLVER_SRC_BOUNDARY_HPP_

#include <vector>

enum class BoundaryType {
  eInherit,
  eExternalTransparent,
  eExternalNoSlipWall,
  eExternalFree
};

class Boundary {
 public:
  Boundary(BoundaryType type_);

  void addPoint(size_t i, size_t j);

 private:
  BoundaryType type;

  std::vector<size_t> i;
  std::vector<size_t> j;
};

#endif
