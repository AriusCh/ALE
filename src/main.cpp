#include <Eigen/Core>
#include <format>
#include <iostream>
#include <memory>
#include <string>

#include "logger.hpp"
#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

#define PROBLEM triplePointShock

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  extern std::shared_ptr<Problem> PROBLEM;
  std::shared_ptr<FEMALEMethod> mtd =
      std::make_shared<FEMALEMethod>("2-48", PROBLEM, 112, 48, 2);
  Simulation sim(mtd);
  sim.run();

  return 0;
}
