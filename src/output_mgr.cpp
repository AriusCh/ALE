#include "output_mgr.hpp"

Writer::Writer() {
  if (!std::filesystem::is_directory("output")) {
    std::filesystem::create_directory("output");
    logger.Log("Created directory \"output\"", LogLevel::eInfo);
  }
  outputDirPath = std::filesystem::current_path() / "output";
}

void Writer::dumpData(std::shared_ptr<GridALE> grid) {}
