#include <Eigen/Core>
#include <memory>

#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  const Problem &problem = Problems::sodTest;
  // const Problem problem = Problem::createRiemannProblem1Dx(
  //     "sod-test-strong", 0.0, 1.0, 0.05, std::deque<double>{}, 0.0,
  //     0.125, 7.0, 0.0, 0.125, 0.1, 0.5, 5.0 / 3.0);
  std::unique_ptr<FEMALEMethod> mtd =
      std::make_unique<FEMALEMethod>(problem, 100, 1, 6, false);
  Simulation sim(std::move(mtd));
  sim.run();

  return 0;
}
