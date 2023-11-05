#include "simulation.hpp"

#include <cmath>
#include <format>

Simulation::Simulation(std::shared_ptr<Method> method,
                       std::shared_ptr<Problem> problem)
    : method(std::move(method)), problem(std::move(problem)) {}
void Simulation::run() {
  t = method->tmin;

  int iterationNum = 0;
  logSimulationStart();
  auto simStart = std::chrono::high_resolution_clock::now();
  while (t < method->tmax) {
    double dt = method->calcdt();

    if (t + dt > method->tmax) {
      dt = method->tmax - t;
    }

    auto start = std::chrono::high_resolution_clock::now();

    method->calc(dt);

    double calcTime = std::chrono::duration<double>(
                          std::chrono::high_resolution_clock::now() - start)
                          .count();

    logSuccessfulIteration(iterationNum, t, dt, calcTime);

    t += dt;
    iterationNum++;
  }
  double simulationTime =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -
                                    simStart)
          .count();
  logSimulationEnd(simulationTime);
}

void Simulation::logSuccessfulIteration(int iterationNumber, double t,
                                        double dt, double calcTime) const {
  std::string message = std::format(
      "PROBLEM: {} ITERATION: {:6} t: {:6.6e} dt: {:6.6e} STEP TIME: {:.6f}",
      problem->name, iterationNumber, t, dt, calcTime);

  logger.Log(message, LogLevel::eGeneral);
}
void Simulation::logSimulationStart() const {
  std::string message =
      std::format("PROBLEM: {} STARTING SIMULATION", problem->name);

  logger.Log(message, LogLevel::eGeneral);
}
void Simulation::logSimulationEnd(double simTime) const {
  std::string message = std::format(
      "PROBLEM: {} SIMULATION COMPLETE. TIME: {:.6f}", problem->name, simTime);

  logger.Log(message, LogLevel::eGeneral);
}
