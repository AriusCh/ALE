#ifndef ALE_SOLVER_SRC_METHOD_HPP_
#define ALE_SOLVER_SRC_METHOD_HPP_

#include "grid.hpp"
#include "problem.hpp"

enum class MethodType { eALE };

class Method {
 public:
  Method(MethodType type);
  Method(Method const &rhs) = delete;
  Method(Method &&rhs) = default;

  Method &operator=(Method const &rhs) = delete;
  Method &operator=(Method &&rhs) = default;

  ~Method();

 public:
  virtual void calc(double dt) = 0;
  virtual double calcdt() const = 0;

  MethodType getType() const;

 protected:
  double tmin;
  double tmax;

 private:
  MethodType type;
};

class MethodALE : public Method {
 public:
  MethodALE(std::unique_ptr<Problem> problem);

 public:
  virtual void calc(double dt);
  virtual double calcdt() const;

 private:
  void initializeMethod();

  void resolveBoundaries();

 private:
  std::vector<std::unique_ptr<GridALE>> grids;
  std::vector<std::unique_ptr<Boundary>> boundaries;
};

#endif
