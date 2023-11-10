#include <iostream>
#include <memory>

#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  extern std::shared_ptr<Problem> blastWave;
  std::shared_ptr<MethodALE> method(
      std::make_shared<MethodALE>(blastWave, 500, 500, 0.0, 0.0, 0.01, 2));
  std::shared_ptr<Simulation> sim(
      std::make_shared<Simulation>(method, blastWave));
  sim->run();

  return 0;
}
