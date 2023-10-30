#include <iostream>
#include <memory>

#include "method.hpp"
#include "problem.hpp"

int main(int argc, char *argv[]) {
  std::unique_ptr<Problem> problem = std::make_unique<RiemannProblem1Dx>(
      0.0, 1.0, 0.2, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 0.5, 1.4);

  MethodALE method(std::move(problem));

  return 0;
}
