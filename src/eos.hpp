#ifndef ALE_SOLVER_SRC_EOS_HPP_
#define ALE_SOLVER_SRC_EOS_HPP_

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
  EOSIdeal(double gamma_);
  EOSIdeal(EOSIdeal const &rhs) = default;
  EOSIdeal(EOSIdeal &&rhs) = default;

  EOSIdeal &operator=(EOSIdeal const &rhs) = default;
  EOSIdeal &operator=(EOSIdeal &&rhs) = default;
  virtual ~EOSIdeal() = default;

 public:
  virtual double getp(double rho, double e) const override;
  virtual double gete(double rho, double p) const override;
  virtual double getc(double rho, double p) const override;
  virtual double gets(double rho, double p) const override;

  double getGamma() const;

 private:
  double gamma;
};

#endif
