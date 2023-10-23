#include "logger.hpp"

#include <iostream>

Logger::Logger() {
#ifdef NDEBUG
  minLogLevel = LogLevel::eWarning;
#else
  minLogLevel = LogLevel::eInfo;
#endif
}
Logger::Logger(LogLevel minLogLevel_) : minLogLevel(minLogLevel_) {}

void Logger::Log(std::string message, LogLevel logLevel) const {
  if (logLevel < minLogLevel) {
    return;
  }

  std::cout << "[";
  switch (logLevel) {
    case LogLevel::eInfo:
      std::cout << "INFO";
    case LogLevel::eWarning:
      std::cout << "WARNING";
    case LogLevel::eError:
      std::cout << "ERROR";
  }
  std::cout << "] " << message << std::endl;
}

LogLevel Logger::getMinLogLevel() const { return minLogLevel; }
void Logger::setMinLogLevel(LogLevel minLogLevel_) {
  minLogLevel = minLogLevel_;
}
