#include <Eigen/Core>
#include <format>
#include <iostream>
#include <memory>
#include <string>

#include "logger.hpp"
#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  const Problem &problem = Problems::sodTest;
  std::unique_ptr<FEMALEMethod> mtd =
      std::make_unique<FEMALEMethod>("1-2000", problem, 2000, 1, 1);
  Simulation sim(std::move(mtd));
  sim.run();

  return 0;
}
