#include "logger.hpp"

#include <iostream>
#include <syncstream>

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

  std::osyncstream ostream(std::cout);
  ostream << "[";
  switch (logLevel) {
    case LogLevel::eInfo:
      ostream << "INFO";
    case LogLevel::eWarning:
      ostream << "WARNING";
    case LogLevel::eError:
      ostream << "ERROR";
  }
  ostream << "] " << message << std::endl;
}

LogLevel Logger::getMinLogLevel() const { return minLogLevel; }
void Logger::setMinLogLevel(LogLevel minLogLevel_) {
  minLogLevel = minLogLevel_;
}
