#include <iostream>
#include <memory>

#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  extern std::shared_ptr<Problem> sodTest;

  std::shared_ptr<MethodALE> method(
      std::make_shared<MethodALE>(sodTest, 100, 100, 1.0, 1));

  std::shared_ptr<Simulation> sim(std::make_shared<Simulation>(method));

  return 0;
}
