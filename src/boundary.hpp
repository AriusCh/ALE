#ifndef ALE_SOLVER_SRC_BOUNDARY_HPP_
#define ALE_SOLVER_SRC_BOUNDARY_HPP_

#include <vector>

enum class BoundaryType {
  eExternalTransparent,
  eExternalNoSlipWall,
  eExternalFree
};

class Boundary {
 public:
  Boundary(BoundaryType type_);

 public:
  BoundaryType getType() const;

 private:
  BoundaryType type;
};

#endif
