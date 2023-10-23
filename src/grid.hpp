#ifndef ALE_SOLVER_SRC_GRID_HPP_
#define ALE_SOLVER_SRC_GRID_HPP_

#include <memory>
#include <vector>

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
  GridALE(size_t sizeX_, size_t sizeY_,
          std::unique_ptr<PrimitiveRectangle> rect);
  GridALE(GridALE const &rhs) = default;
  GridALE(GridALE &&rhs) = default;

  GridALE &operator=(GridALE const &rhs) = default;
  GridALE &operator=(GridALE &&rhs) = default;

  ~GridALE() = default;

 public:
  void calc() override;

 private:
  std::vector<std::vector<double>> x;
  std::vector<std::vector<double>> y;
  std::vector<std::vector<double>> rho;
  std::vector<std::vector<double>> u;
  std::vector<std::vector<double>> v;
  std::vector<std::vector<double>> p;
  std::vector<std::vector<double>> m;

  std::unique_ptr<Boundary> leftBcs;
  std::unique_ptr<Boundary> rightBcs;
  std::unique_ptr<Boundary> bottomBcs;
  std::unique_ptr<Boundary> topBcs;
};

#endif
