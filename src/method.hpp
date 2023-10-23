#ifndef ALE_SOLVER_METHOD_HPP_
#define ALE_SOLVER_METHOD_HPP_

#include "grid.hpp"

enum class MethodType { eALE };

class Method {
 public:
  Method();
  Method(Method const &rhs);
  Method(Method &&rhs);

  Method &operator=(Method const &rhs);
  Method &operator=(Method &&rhs);

  ~Method();

 private:
  MethodType type;
};

class MethodALE : public Method {
 public:
  MethodALE();

  private:
  GridALE grid;
};

#endif
