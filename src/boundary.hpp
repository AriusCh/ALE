#ifndef ALE_SOLVER_SRC_BOUNDARY_HPP_
#define ALE_SOLVER_SRC_BOUNDARY_HPP_

#include <vector>

enum class BoundaryType {
  eExternalFree = 0,
  eExternalTransparent = 1,
  eExternalNoSlipWall = 2,
};

#endif
