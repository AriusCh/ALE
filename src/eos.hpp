#ifndef ALE_SOLVER_SRC_EOS_HPP_
#define ALE_SOLVER_SRC_EOS_HPP_

#include <cmath>
class EOS {
 public:
  virtual ~EOS() = 0;

 public:
  virtual double getp(double rho, double e) const = 0;
  virtual double gete(double rho, double p) const = 0;
  virtual double getc(double rho, double p) const = 0;
  virtual double gets(double rho, double p) const = 0;
};

class EOSIdeal : public EOS {
 public:
  EOSIdeal(double gamma_) : gamma(gamma_) {}
  EOSIdeal(EOSIdeal const &rhs) = default;
  EOSIdeal(EOSIdeal &&rhs) = default;

  EOSIdeal &operator=(EOSIdeal const &rhs) = default;
  EOSIdeal &operator=(EOSIdeal &&rhs) = default;
  virtual ~EOSIdeal() = default;

 public:
  inline virtual double getp(double rho, double e) const override {
    return (gamma - 1.) * rho * e;
  }
  inline virtual double gete(double rho, double p) const override {
    return p / (gamma - 1.) / rho;
  }
  inline virtual double getc(double rho, double p) const override {
    return std::sqrt(gamma * p / rho);
  }
  inline virtual double gets(double rho, double p) const override {
    return p / std::pow(rho, gamma);
  }

  inline double getGamma() const { return gamma; };

 private:
  double gamma;
};

#endif
