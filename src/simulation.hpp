#ifndef ALE_SOLVER_SRC_SIMULATION_HPP_
#define ALE_SOLVER_SRC_SIMULATION_HPP_

#include <memory>
#include <queue>

#include "method.hpp"

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

  // Step time relaxation
  static constexpr size_t stepTimesMaxSize = 10;
  std::queue<double> stepTimes;
  double stepTimesSum = 0.0;

  Logger logger;
};

#endif
