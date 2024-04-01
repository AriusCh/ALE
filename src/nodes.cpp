#include "nodes.hpp"

#include <cmath>

double lobattoBasis1D(double x, size_t order, size_t k) {
  assert(order <= lobattoOrderMax);
  assert(x >= 0.0 && x <= 1.0);
  assert(k < order + 1);
  if (order == 0) {
    return 1.0;
  }
  double numerator = 1.0;
  double denominator = 1.0;
  const size_t indMin = getLobattoStartIndex(order);
  const size_t indJ = indMin + k;
  const size_t indMax = indMin + order + 1;

  for (size_t i = indMin; i < indJ; i++) {
    numerator *= x - lobattoAbscissas[i];
    assert(lobattoAbscissas[indJ] != lobattoAbscissas[i]);
    denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[i];
  }
  for (size_t i = indJ + 1; i < indMax; i++) {
    numerator *= x - lobattoAbscissas[i];
    assert(lobattoAbscissas[indJ] != lobattoAbscissas[i]);
    denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[i];
  }
  return numerator / denominator;
}

double lobattoBasis1Ddx(double x, size_t order, size_t k) {
  assert(order <= lobattoOrderMax && order > 0);
  assert(x >= 0.0 && x <= 1.0);
  assert(k < order + 1);
  double output = 0.0;
  const size_t indMin = (order - 1) * (order + 2) / 2;
  const size_t indJ = indMin + k;
  const size_t indMax = indMin + order + 1;
  for (size_t i = indMin; i < indJ; i++) {
    double numerator = 1.0;
    double denominator = 1.0;
    for (size_t j = indMin; j < i; j++) {
      numerator *= x - lobattoAbscissas[j];
      assert(lobattoAbscissas[indJ] != lobattoAbscissas[j]);
      denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[j];
    }
    for (size_t j = i + 1; j < indJ; j++) {
      numerator *= x - lobattoAbscissas[j];
      assert(lobattoAbscissas[indJ] != lobattoAbscissas[j]);
      denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[j];
    }
    for (size_t j = indJ + 1; j < indMax; j++) {
      numerator *= x - lobattoAbscissas[j];
      assert(lobattoAbscissas[indJ] != lobattoAbscissas[j]);
      denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[j];
    }
    assert(lobattoAbscissas[indJ] != lobattoAbscissas[i]);
    output += numerator /
              (denominator * (lobattoAbscissas[indJ] - lobattoAbscissas[i]));
  }
  for (size_t i = indJ + 1; i < indMax; i++) {
    double numerator = 1.0;
    double denominator = 1.0;
    for (size_t j = indMin; j < indJ; j++) {
      numerator *= x - lobattoAbscissas[j];
      assert(lobattoAbscissas[indJ] != lobattoAbscissas[j]);
      denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[j];
    }
    for (size_t j = indJ + 1; j < i; j++) {
      numerator *= x - lobattoAbscissas[j];
      assert(lobattoAbscissas[indJ] != lobattoAbscissas[j]);
      denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[j];
    }
    for (size_t j = i + 1; j < indMax; j++) {
      numerator *= x - lobattoAbscissas[j];
      assert(lobattoAbscissas[indJ] != lobattoAbscissas[j]);
      denominator *= lobattoAbscissas[indJ] - lobattoAbscissas[j];
    }
    assert(lobattoAbscissas[indJ] != lobattoAbscissas[i]);
    output += numerator /
              (denominator * (lobattoAbscissas[indJ] - lobattoAbscissas[i]));
  }
  return output;
}

double legendreBasis1D(double x, size_t order, size_t k) {
  assert(order + 1 <= legendreOrderMax);
  assert(x >= 0.0 && x <= 1.0);
  assert(k < order + 1);

  double numerator = 1.0;
  double denominator = 1.0;
  const size_t indMin = getLegendreStartIndex(order + 1);
  const size_t indJ = indMin + k;
  const size_t indMax = indMin + order + 1;

  for (size_t i = indMin; i < indJ; i++) {
    numerator *= x - legendreAbscissas[i];
    assert(legendreAbscissas[indJ] != legendreAbscissas[i]);
    denominator *= legendreAbscissas[indJ] - legendreAbscissas[i];
  }
  for (size_t i = indJ + 1; i < indMax; i++) {
    numerator *= x - legendreAbscissas[i];
    assert(legendreAbscissas[indJ] != legendreAbscissas[i]);
    denominator *= legendreAbscissas[indJ] - legendreAbscissas[i];
  }

  return numerator / denominator;
}

double bernsteinBasis1D(double x, size_t order, size_t k) {
  assert(k <= order);
  if (order == 0) {
    return 1.0;
  }
  double kmin = std::min(k, order - k);
  double numerator = 1.0;
  double denominator = 1.0;
  for (size_t i = 1; i <= kmin; i++) {
    numerator *= order - kmin + i;
    denominator *= i;
  }
  double coeff = numerator / denominator;
  return coeff * std::pow(x, k) * std::pow(1.0 - x, order - k);
}
double bernsteinBasis1Ddx(double x, size_t order, size_t k) {
  assert(k <= order);
  if (order == 0) {
    return 0.0;
  }
  double kmin = std::min(k, order - k);
  double numerator = 1.0;
  double denominator = 1.0;
  for (size_t i = 1; i <= kmin; i++) {
    numerator *= order - kmin + i;
    denominator *= i;
  }
  double coeff = numerator / denominator;
  double leftPart = 0.0;
  double rightPart = 0.0;
  if (k > 0) {
    leftPart = k * std::pow(x, k - 1) * std::pow(1.0 - x, order - k);
  }
  if (k < order) {
    rightPart = (order - k) * std::pow(x, k) * std::pow(1.0 - x, order - k - 1);
  }
  return coeff * (leftPart - rightPart);
}
