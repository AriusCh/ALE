#include "simulation.hpp"

#include <cmath>
#include <format>

Simulation::Simulation(std::shared_ptr<Method> method)
    : method(std::move(method)) {}
void Simulation::run() {
  method->t = method->tmin;

  int iterationNum = 0;
  logSimulationStart();
  method->dumpGrid();
  auto simStart = std::chrono::high_resolution_clock::now();
  while (method->t < method->tmax) {
    method->calcdt();
    if (method->t + method->dt > method->tmax) {
      method->dt = method->tmax - method->t;
    }

    auto start = std::chrono::high_resolution_clock::now();

    method->calc();

    double calcTime = std::chrono::duration<double>(
                          std::chrono::high_resolution_clock::now() - start)
                          .count();

    logSuccessfulIteration(iterationNum, method->t, method->dt, calcTime);

    iterationNum++;
  }
  double simulationTime =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -
                                    simStart)
          .count();
  logSimulationEnd(simulationTime);
  method->dumpGrid();
}

void Simulation::logSuccessfulIteration(int iterationNumber, double t,
                                        double dt, double calcTime) const {
  std::string message = std::format(
      "PROBLEM: {} ITERATION: {:6} t: {:6.6e} dt: {:6.6e} STEP TIME: {:.6f}",
      method->problemName, iterationNumber, t, dt, calcTime);

  logger.log(message, LogLevel::eGeneral);
}
void Simulation::logSimulationStart() const {
  std::string message =
      std::format("PROBLEM: {} STARTING SIMULATION", method->problemName);

  logger.log(message, LogLevel::eGeneral);
}
void Simulation::logSimulationEnd(double simTime) const {
  std::string message =
      std::format("PROBLEM: {} SIMULATION COMPLETE. TIME: {:.6f}",
                  method->problemName, simTime);

  logger.log(message, LogLevel::eGeneral);
}
