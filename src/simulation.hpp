#ifndef ALE_SOLVER_SRC_SIMULATION_HPP_
#define ALE_SOLVER_SRC_SIMULATION_HPP_

class Simulation {
 public:
  Simulation();
  Simulation(Simulation const &rhs);
  Simulation(Simulation &&rhs);

  Simulation &operator=(Simulation const &rhs);
  Simulation &operator=(Simulation &&rhs);

  ~Simulation();

};

#endif
