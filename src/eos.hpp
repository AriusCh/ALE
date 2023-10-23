#ifndef ALE_SOLVER_SRC_EOS_HPP_
#define ALE_SOLVER_SRC_EOS_HPP_

class EOS {
 public:
  virtual double getP() = 0;
  virtual double getE() = 0;
  virtual double getC() = 0;
};

#endif

