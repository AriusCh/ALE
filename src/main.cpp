#include <Eigen/Core>
#include <format>
#include <iostream>
#include <memory>
#include <string>

#include "logger.hpp"
#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

#define PROBLEM blastWave2D

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  extern std::shared_ptr<Problem> PROBLEM;
  std::shared_ptr<FEMALEMethod> mtd =
      std::make_shared<FEMALEMethod>(PROBLEM->name, PROBLEM, 50, 50, 1);
  Simulation sim(mtd);
  sim.run();

  return 0;
}
