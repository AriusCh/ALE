#ifndef ALE_SOLVER_SRC_GRID_HPP_
#define ALE_SOLVER_SRC_GRID_HPP_

#include <memory>
#include <vector>

#include "boundary.hpp"
#include "eos.hpp"
#include "logger.hpp"
#include "primitive.hpp"

enum class GridType { eGridALE };

class Grid {
 public:
  Grid(GridType type_, size_t sizeX_, size_t sizeY_);
  Grid(Grid const &rhs) = default;
  Grid(Grid &&rhs) = default;

  Grid &operator=(Grid const &rhs) = default;
  Grid &operator=(Grid &&rhs) = default;

  ~Grid() = default;

  size_t getSizeX() const;
  size_t getSizeY() const;

  virtual void calc() = 0;

 public:
 protected:
  size_t sizeX;
  size_t sizeY;

 private:
  GridType type;
  Logger logger;
};

class GridALE : public Grid {
 public:
  GridALE(const std::vector<std::vector<double>> &x,
          const std::vector<std::vector<double>> &y,
          const std::vector<std::vector<double>> &rho,
          const std::vector<std::vector<double>> &p,
          const std::vector<std::vector<double>> &u,
          const std::vector<std::vector<double>> &v,
          std::unique_ptr<EOS> &&eos);
  GridALE(GridALE const &rhs) = delete;
  GridALE(GridALE &&rhs) = default;

  GridALE &operator=(GridALE const &rhs) = delete;
  GridALE &operator=(GridALE &&rhs) = default;

  ~GridALE() = default;

 public:
  void calc() override;

 private:
  std::vector<std::vector<double>> x;    // Node x coordinates
  std::vector<std::vector<double>> y;    // Node y coordinates
  std::vector<std::vector<double>> rho;  // Cell density
  std::vector<std::vector<double>> p;    // Cell pressure
  std::vector<std::vector<double>> u;    // Node velocity in Ox direction
  std::vector<std::vector<double>> v;    // Node velocity in Oy direction
  std::vector<std::vector<double>> m;    // Node mass
  std::vector<std::vector<double>> E;    // Cell total specific energy

  std::shared_ptr<Boundary> leftBcs;
  std::shared_ptr<Boundary> rightBcs;
  std::shared_ptr<Boundary> bottomBcs;
  std::shared_ptr<Boundary> topBcs;

  std::unique_ptr<EOS> eos;
};

#endif
