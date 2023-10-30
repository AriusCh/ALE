#include "simulation.hpp"

Simulation::Simulation(std::unique_ptr<Method> &&method)
    : method(std::move(method)) {}
void Simulation::run() {
    while 
}
