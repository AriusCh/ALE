#ifndef ALE_SOLVER_SRC_OUTPUT_MGR_HPP_
#define ALE_SOLVER_SRC_OUTPUT_MGR_HPP_

#include <filesystem>

#include "grid.hpp"
#include "logger.hpp"

class Writer {
 public:
  Writer();

  void dumpData(std::shared_ptr<GridALE> grid);

 private:
  std::filesystem::path outputDirPath;

  Logger logger;
};

#endif
