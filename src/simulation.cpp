#include "simulation.hpp"

#include <cmath>
#include <format>

Simulation::Simulation(std::unique_ptr<Method> method)
    : method(std::move(method)) {}
void Simulation::run() {
  method->t = method->problem.tmin;
  auto tOut = method->problem.tOut;

  int iterationNum = 0;
  logSimulationStart();
  method->dumpData();
  method->dumpGrid();
  auto simStart = std::chrono::high_resolution_clock::now();
  while (method->t < method->problem.tmax) {
    method->calcdt();
    double tmpDt = method->dt;
    if (!tOut.empty()) {
      double tNext = tOut.front();
      if (method->t + method->dt > tNext) {
        method->dt = tNext - method->t;
      }
    }
    if (method->t + method->dt > method->problem.tmax) {
      method->dt = method->problem.tmax - method->t;
    }

    auto start = std::chrono::high_resolution_clock::now();

    method->calc();

    double calcTime = std::chrono::duration<double>(
                          std::chrono::high_resolution_clock::now() - start)
                          .count();

    double remTime = (method->problem.tmax - method->t) / method->dt * calcTime;

    logSuccessfulIteration(iterationNum, method->t, method->dt, calcTime,
                           remTime);

    if (!tOut.empty()) {
      double tNext = tOut.front();
      if (method->t >= tNext) {
        method->dumpData();
        method->dumpGrid();
        tOut.pop_front();
        method->dt = tmpDt;
      }
    }
    iterationNum++;
  }
  double simulationTime =
      std::chrono::duration<double>(std::chrono::high_resolution_clock::now() -
                                    simStart)
          .count();
  logSimulationEnd(simulationTime);
  method->dumpData();
  method->dumpGrid();
}

void Simulation::logSuccessfulIteration(int iterationNumber, double t,
                                        double dt, double calcTime,
                                        double remTime) const {
  std::string message = std::format(
      "PROBLEM: {} ITERATION: {:6} t: {:6.6e} dt: {:6.6e} STEP TIME: {:.6f} "
      "REMAINING TIME: {:.0f}",
      method->problem.name, iterationNumber, t, dt, calcTime, remTime);

  logger.log(message, LogLevel::eGeneral);
}
void Simulation::logSimulationStart() const {
  std::string message =
      std::format("PROBLEM: {} STARTING SIMULATION", method->problem.name);

  logger.log(message, LogLevel::eGeneral);
}
void Simulation::logSimulationEnd(double simTime) const {
  std::string message =
      std::format("PROBLEM: {} SIMULATION COMPLETE. TIME: {:.6f}",
                  method->problem.name, simTime);

  logger.log(message, LogLevel::eGeneral);
}
