#include "logger.hpp"

#include <iostream>
#include <syncstream>

Logger::Logger() {
#ifdef NDEBUG
  minLogLevel = LogLevel::eGeneral;
#else
  minLogLevel = LogLevel::eInfo;
#endif
}
Logger::Logger(LogLevel minLogLevel_) : minLogLevel(minLogLevel_) {}

void Logger::log(std::string message, LogLevel logLevel) const {
  if (logLevel < minLogLevel) {
    return;
  }

  std::osyncstream ostream(std::cout);
  ostream << "[";
  switch (logLevel) {
    case LogLevel::eInfo:
      ostream << "INFO";
      break;
    case LogLevel::eGeneral:
      ostream << "GENERAL";
      break;
    case LogLevel::eWarning:
      ostream << "WARNING";
      break;
    case LogLevel::eError:
      ostream << "ERROR";
      break;
  }
  ostream << "] " << message << std::endl;
}

LogLevel Logger::getMinLogLevel() const { return minLogLevel; }
void Logger::setMinLogLevel(LogLevel minLogLevel_) {
  minLogLevel = minLogLevel_;
}
