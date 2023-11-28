#include <iostream>
#include <memory>

#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

#define PROBLEM toro3x;

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  extern std::shared_ptr<Problem> PROBLEM;
  auto problem = PROBLEM;
  std::shared_ptr<MethodALE> method(
      std::make_shared<MethodALE>(problem, 100, 100, 0.0, 0.0, 0.1, 1));
  std::shared_ptr<Simulation> sim(
      std::make_shared<Simulation>(method, problem));
  sim->run();

  return 0;
}
