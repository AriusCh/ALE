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

class EOSMGAlPrecise6 : public EOS {
 public:
  EOSMGAlPrecise6()
      : rho0(2750.), p0(560.964e9), a(1.12657), b(0.975511), ph(15.e9) {}

 public:
  double getp(double rho, double e) const override;
  double gete(double rho, double p) const override;
  double getc(double rho, double p) const override;
  double gets(double rho, double p) const override;
  double getdpdrho(double rho, double e) const;
  double getdpde(double rho, double e) const;

 private:
  double G(double x) const;
  double GPrime(double x) const;
  // double Gx1(double x) const;
  // double Gx2(double x) const;
  double pCold(double rho) const;
  double pColdPrime(double rho) const;
  double eCold(double rho) const;
  double eColdPrime(double rho) const;

 private:
  const double rho0, p0, a, b, ph;
};

#endif
