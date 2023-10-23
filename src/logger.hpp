#ifndef ALE_SOLVER_SRC_LOGGER_HPP_
#define ALE_SOLVER_SRC_LOGGER_HPP_

#include <string>

enum class LogLevel { eInfo = 0, eWarning = 1, eError = 2 };

class Logger {
 public:
  Logger();
  Logger(LogLevel minLogLevel_);
  Logger(Logger const &rhs) = default;
  Logger(Logger &&rhs) = default;

  Logger &operator=(Logger const &rhs) = default;
  Logger &operator=(Logger &&rhs) = default;

  ~Logger() = default;

  void Log(std::string message, LogLevel logLevel = LogLevel::eInfo) const;

  LogLevel getMinLogLevel() const;
  void setMinLogLevel(LogLevel minLogLevel_);

 private:
  LogLevel minLogLevel = LogLevel::eInfo;
};

#endif
