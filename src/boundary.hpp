#ifndef ALE_SOLVER_SRC_BOUNDARY_HPP_
#define ALE_SOLVER_SRC_BOUNDARY_HPP_

#include <vector>

enum class BoundaryType {
  eFree = 0,
  eWall = 1,
  eNoSlipWall = 2,
};

#endif
