#ifndef ALE_SOLVER_SRC_EOS_HPP_
#define ALE_SOLVER_SRC_EOS_HPP_

class EOS {
 public:
  virtual double getp(double rho, double e) const = 0;
  virtual double gete(double rho, double p) const = 0;
  virtual double getc(double rho, double p) const = 0;
};

class EOSIdeal : public EOS {
 public:
  EOSIdeal(double gamma_);

 public:
  virtual double getp(double rho, double e) const override;
  virtual double gete(double rho, double p) const override;
  virtual double getc(double rho, double p) const override;

  double getGamma() const;

 private:
  double gamma;
};

#endif
