#include "eos.hpp"

#include <cmath>

EOS::~EOS() {}

EOSIdeal::EOSIdeal(double gamma_) : gamma(gamma_) {}
double EOSIdeal::getp(double rho, double e) const {
  return (gamma - 1.) * rho * e;
}
double EOSIdeal::gete(double rho, double p) const {
  return p / (gamma - 1.) / rho;
}
double EOSIdeal::getc(double rho, double p) const {
  return std::sqrt(gamma * p / rho);
}
double EOSIdeal::gets(double rho, double p) const {
  return p / std::pow(rho, gamma);
}
double EOSIdeal::getGamma() const { return gamma; }
