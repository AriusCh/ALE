#ifndef ALE_SOLVER_SRC_SIMULATION_HPP_
#define ALE_SOLVER_SRC_SIMULATION_HPP_

#include <memory>

#include "method.hpp"
#include "output_mgr.hpp"

class Simulation {
 public:
  Simulation(std::unique_ptr<Method> method);
  Simulation(Simulation const &rhs) = delete;
  Simulation(Simulation &&rhs) = default;

  Simulation &operator=(Simulation const &rhs) = delete;
  Simulation &operator=(Simulation &&rhs) = default;

  ~Simulation() = default;

 public:
  void run();

 private:
  void logSuccessfulIteration(int iterationNumber, double t, double dt,
                              double calcTime, double remTime) const;
  void logSimulationStart() const;
  void logSimulationEnd(double simTime) const;

 private:
  std::unique_ptr<Method> method;

  Logger logger;
};

#endif
