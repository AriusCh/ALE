#ifndef ALE_SOLVER_SRC_OUTPUT_MGR_HPP_
#define ALE_SOLVER_SRC_OUTPUT_MGR_HPP_

#include <filesystem>
#include <fstream>
#include <functional>
// #include <future>
#include <string>

#include "logger.hpp"

class Writer {
 public:
  Writer(const std::string &problemName, const std::string &name);
  ~Writer();

  void dumpData(std::function<void(std::ofstream ofs)>,
                std::string filename) const;

 private:
  std::filesystem::path outputDirPath;

  Logger logger;
};

#endif
