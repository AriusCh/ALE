#ifndef ALE_SOLVER_SRC_SIMULATION_HPP_
#define ALE_SOLVER_SRC_SIMULATION_HPP_

#include <memory>

#include "method.hpp"
class Simulation {
 public:
  Simulation(std::unique_ptr<Method> &&method);
  Simulation(Simulation const &rhs);
  Simulation(Simulation &&rhs);

  Simulation &operator=(Simulation const &rhs);
  Simulation &operator=(Simulation &&rhs);

  ~Simulation();

  void run();

 private:
  std::unique_ptr<Method> method;

  double t;
};

#endif
