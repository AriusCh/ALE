#ifndef ALE_SOLVER_SRC_METHOD_HPP_
#define ALE_SOLVER_SRC_METHOD_HPP_

#include <barrier>
#include <future>

#include "grid.hpp"
#include "problem.hpp"

enum class Status { eExit, eStepStart, eCalcdt };

class Method {
 public:
  Method(double tmin_, double tmax_) : t(tmin_), tmin(tmin_), tmax(tmax_) {}
  Method(Method const &rhs) = default;
  Method(Method &&rhs) = default;

  Method &operator=(Method const &rhs) = delete;
  Method &operator=(Method &&rhs) = delete;

  virtual ~Method() = default;

 public:
  virtual void calc(double dt) = 0;
  virtual double calcdt() const = 0;

  // virtual void dumpGrid(std::shared_ptr<Problem> problem) const = 0;

 public:
  double t;

  const double tmin;
  const double tmax;
};

class FEMALEMethod : public Method {
 public:
  FEMALEMethod(std::shared_ptr<Problem> pr_, size_t xCells_, size_t yCells_,
               size_t order_);
  FEMALEMethod(FEMALEMethod const &rhs) = default;
  FEMALEMethod(FEMALEMethod &&rhs) = default;

  FEMALEMethod &operator=(FEMALEMethod const &rhs) = delete;
  FEMALEMethod &operator=(FEMALEMethod &&rhs) = delete;

  virtual ~FEMALEMethod() = default;

 public:
  virtual void calc(double dt);
  virtual double caltdt() const;

 private:
  size_t xCells;
  size_t yCells;
  size_t order;
};

#endif
