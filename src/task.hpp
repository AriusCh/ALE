#ifndef ALE_SOLVER_SRC_TASK_HPP_
#define ALE_SOLVER_SRC_TASK_HPP_

#include "logger.hpp"

class Task {
 public:
  Task();
  Task(Task const &rhs);
  Task(Task &&rhs);

  Task &operator=(Task const &rhs);
  Task &operator=(Task &&rhs);

  ~Task();

 private:
  Logger logger;
};

#endif
