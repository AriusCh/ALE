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
  const Problem &problem = Problems::task5v11;
  {
    std::unique_ptr<FEMALEMethod> mtd =
        std::make_unique<FEMALEMethod>("1-490", problem, 490, 175, 1);
    Simulation sim(std::move(mtd));
    sim.run();
  }
  {
    std::unique_ptr<FEMALEMethod> mtd =
        std::make_unique<FEMALEMethod>("1-700", problem, 700, 250, 1);
    Simulation sim(std::move(mtd));
    sim.run();
  }

  return 0;
}
