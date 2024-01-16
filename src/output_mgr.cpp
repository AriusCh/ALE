#include "output_mgr.hpp"

#include <format>

Writer::Writer(const std::string& problemName) {
  if (!std::filesystem::is_directory("output")) {
    std::filesystem::create_directory("output");
    logger.log("Created directory \"output\"", LogLevel::eInfo);
  }
  if (!std::filesystem::is_directory("output/" + problemName)) {
    std::filesystem::create_directory("output/" + problemName);
    logger.log("Created directory \"output/" + problemName + "\"",
               LogLevel::eInfo);
  }
  outputDirPath = std::filesystem::current_path() / "output" / problemName;
}
Writer::~Writer() {}

void Writer::dumpData(std::function<void(std::ofstream ofs)> func,
                      std::string filename) const {
  std::ofstream ofs(outputDirPath / filename);
  if (!ofs) {
    std::string message = std::format("FAILED TO OPEN FILE: {}",
                                      (outputDirPath / filename).string());
    logger.log(message, LogLevel::eError);
    return;
  }
  std::string message =
      std::format("WRITING TO FILE: {}", (outputDirPath / filename).string());
  logger.log(message, LogLevel::eGeneral);
  func(std::move(ofs));
}
