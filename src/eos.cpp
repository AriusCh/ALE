#include "eos.hpp"

#include <cmath>

EOS::~EOS() {}

double EOSMGAlPrecise6::getp(double rho, double e) const {
  double x = rho / rho0;
  return pCold(rho) + G(x) * rho * (e - eCold(rho));
}

double EOSMGAlPrecise6::gete(double rho, double p) const {
  double x = rho / rho0;
  double _ec = eCold(rho);
  double _pc = pCold(rho);
  return _ec + 1. / rho / G(x) * (p - _pc);
}

double EOSMGAlPrecise6::getc(double rho, double p) const {
  double x = rho / rho0;
  double _G = G(x);
  double G_x = GPrime(x);
  double pc = pCold(rho) / 1.e9;
  double pcx = pColdPrime(rho) / 1.e9;
  double c2 = 1.e9 / rho0 *
              ((p / 1.e9 - pc) / (x * G(x)) * (_G + _G * _G + x * G_x) + pcx);

  return sqrt(c2);
}

double EOSMGAlPrecise6::gets(double rho, double p) const {
  const double x1 = 1., p1 = 0., cv = 1.;
  double x = rho / rho0;
  return cv * log((p - pCold(x)) / (p1 - pCold(1.)) * pow(x1 / x, G(x) + 1.));
}

double EOSMGAlPrecise6::getdpdrho(double rho, double e) const {
  double x = rho / rho0;
  double _G = G(x);
  double _pc = pCold(rho);
  double _pcx = pColdPrime(rho);
  double _ec = eCold(rho);
  // double _ecx = eColdPrime(rho);
  // return _pcx/rho0 + _G*(e-_ec) - rho*_G*_ecx/rho0;
  return 1. / rho0 * (_pcx + _G * rho0 * (e - _ec) - _G / x * _pc);
}

double EOSMGAlPrecise6::getdpde(double rho, [[maybe_unused]] double e) const {
  double _G = G(rho / rho0);
  return rho * _G;
}

double EOSMGAlPrecise6::G([[maybe_unused]] double x) const { return 1.2; }

double EOSMGAlPrecise6::GPrime([[maybe_unused]] double x) const { return 0.; }

double EOSMGAlPrecise6::pCold(double rho) const {
  double x = rho / rho0;
  double pc = 0.;
  if (x >= 1.)
    pc = p0 * x * (pow(x, a) - pow(x, b));
  else {
    double n = p0 * (a - b) / ph;
    pc = ph * (pow(x, n) - 1.);
  }
  return pc;
}

double EOSMGAlPrecise6::pColdPrime(double rho) const {
  double x = rho / rho0;
  double pcx = 0.;
  if (x >= 1.)
    pcx = p0 * ((a + 1.) * pow(x, a) - (b + 1.) * pow(x, b));
  else {
    double n = p0 * (a - b) / ph;
    pcx = p0 * (a - b) * pow(x, n - 1.);
  }
  return pcx;
}

double EOSMGAlPrecise6::eCold(double rho) const {
  double x = rho / rho0;
  double ec = 0.;
  if (x >= 1.)
    ec = 2.03987e8 * (pow(x, a) / a - pow(x, b) / b);
  else {
    double n = p0 * (a - b) / ph;
    ec = 2.03987e8 * (a - b) *
         (1. / n * (pow(x, n - 1.) / (n - 1.) + 1. / x) - 1. / (n - 1.) -
          1. / a / b);
  }
  return ec;
}

double EOSMGAlPrecise6::eColdPrime(double rho) const {
  double x = rho / rho0;
  double ecx = 0.;
  if (x >= 1.)
    ecx = 2.03987e8 * (pow(x, a - 1.) - pow(x, b - 1.));
  else {
    double n = p0 * (a - b) / ph;
    ecx = 2.03987e17 / p0 * (a - b) * (1. / n * (pow(x, n - 2.) - 1. / x / x));
  }
  return ecx;
}
