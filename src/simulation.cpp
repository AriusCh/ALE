#include "simulation.hpp"

#include <cmath>
#include <format>

Simulation::Simulation(std::shared_ptr<Method> method)
    : method(std::move(method)) {}
void Simulation::run() {
  t = method->tmin;

  int iterationNum = 0;
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
}

void Simulation::logSuccessfulIteration(int iterationNumber, double t,
                                        double dt, double calcTime) const {
  std::string message =
      std::format("ITERATION: {:6} t: {:6.6e} dt: {:6.6e} STEP TIME: {:.6f}",
                  iterationNumber, t, dt, calcTime);

  logger.Log(message, LogLevel::eInfo);
}
