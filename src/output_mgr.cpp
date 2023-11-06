#include "output_mgr.hpp"

#include <format>

Writer::Writer(std::string problemName) {
  if (!std::filesystem::is_directory("output")) {
    std::filesystem::create_directory("output");
    logger.Log("Created directory \"output\"", LogLevel::eInfo);
  }
  if (!std::filesystem::is_directory("output/" + problemName)) {
    std::filesystem::create_directory("output/" + problemName);
    logger.Log("Created directory \"output/" + problemName + "\"",
               LogLevel::eInfo);
  }
  outputDirPath = std::filesystem::current_path() / "output" / problemName;
}
Writer::~Writer() {
  for (auto &fut : futures) {
    fut.get();
  }
}

void Writer::dumpData(std::function<void(std::ofstream ofs)> func,
                      std::string filename) const {
  std::ofstream ofs(outputDirPath / filename);
  if (!ofs) {
    std::string message = std::format("FAILED TO OPEN FILE: {}",
                                      (outputDirPath / filename).string());
    logger.Log(message, LogLevel::eError);
    return;
  }
  std::string message =
      std::format("WRITING TO FILE: {}", (outputDirPath / filename).string());
  logger.Log(message, LogLevel::eInfo);
  std::future<void> fut = std::async(func, std::move(ofs));
  fut.wait();
  futures.push_back(std::move(fut));
}
