#ifndef ALE_SOLVER_SRC_GRID_MANAGER_HPP_
#define ALE_SOLVER_SRC_GRID_MANAGER_HPP_

#include <memory>
#include <vector>

#include "grid.hpp"

class GridManager {
  GridManager(std::vector<std::unique_ptr<Grid>> &&grids);

 private:
  std::vector<std::unique_ptr<Grid>> grids;
};

#endif
