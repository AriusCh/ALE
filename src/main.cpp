#include <Eigen/Core>
#include <memory>

#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  const Problem &problem = Problems::sodTest;
  std::unique_ptr<FEMALEMethod> mtd =
      std::make_unique<FEMALEMethod>(problem, 100, 1, 1, false);
  Simulation sim(std::move(mtd));
  sim.run();

  return 0;
}
