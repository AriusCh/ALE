#ifndef ALE_SOLVER_SRC_GRID_HPP_
#define ALE_SOLVER_SRC_GRID_HPP_

#include <functional>
#include <memory>
#include <vector>

#include "boundary.hpp"
#include "eos.hpp"
#include "logger.hpp"

enum class GridType { eFEMALEGrid };

// class GridALE {
//  public:
//   GridALE(int sizeX, int sizeY, double xmin, double xmax, double ymin,
//           double ymax, std::function<double(double, double)> uInit,
//           std::function<double(double, double)> vInit,
//           std::function<double(double, double)> rhoInit,
//           std::function<double(double, double)> pInit,
//           std::shared_ptr<EOS> eos);
//   GridALE(GridALE const &rhs) = default;
//   GridALE(GridALE &&rhs) = default;
//
//   GridALE &operator=(GridALE const &rhs) = delete;
//   GridALE &operator=(GridALE &&rhs) = delete;
//
//   ~GridALE() = default;
//
//  public:
//   /* Calculate volume of cell i, j of a grid (x, y) */
//   static double getV(int i, int j, std::vector<std::vector<double>> &x,
//                      std::vector<std::vector<double>> &y);
//
//  protected:
//   void populateCellEnergy();
//   void populateNodeMass();
//
//  public:
//   const int sizeX, sizeY;
//
//   std::vector<std::vector<double>> x;    // Node x coordinates
//   std::vector<std::vector<double>> y;    // Node y coordinates
//   std::vector<std::vector<double>> u;    // Node velocity in Ox direction
//   std::vector<std::vector<double>> v;    // Node velocity in Oy direction
//   std::vector<std::vector<double>> m;    // Node mass
//   std::vector<std::vector<double>> rho;  // Cell density
//   std::vector<std::vector<double>> p;    // Cell pressure
//   std::vector<std::vector<double>> E;    // Cell total specific energy
//
//   std::shared_ptr<EOS> eos;
// };

// class Grid {
//  public:
//   Grid(GridType type_) : type(type_) {}
//   Grid(Grid const &rhs) = default;
//   Grid(Grid &&rhs) = default;
//
//   Grid &operator=(Grid const &rhs) = default;
//   Grid &operator=(Grid &&rhs) = default;
//
//   virtual ~Grid() = 0;
//
//  public:
//   GridType type;
// };
//
// class FEMALEGrid : Grid {
//   FEMALEGrid();
//   FEMALEGrid(FEMALEGrid const &rhs) = default;
//   FEMALEGrid(FEMALEGrid &&rhs) = default;
//
//   FEMALEGrid &operator=(FEMALEGrid const &rhs) = default;
//   FEMALEGrid &operator=(FEMALEGrid &&rhs) = default;
//
//   virtual ~FEMALEGrid() = 0;
// };

#endif
