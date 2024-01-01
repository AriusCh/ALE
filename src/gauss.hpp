#ifndef ALE_SOLVER_SRC_GAUSS_HPP_
#define ALE_SOLVER_SRC_GAUSS_HPP_

#include <array>
#include <cassert>
#include <functional>

// Legendre nodes
constexpr std::array<double, 1> legendreX1{0.0};
constexpr std::array<double, 2> legendreX2{-0.5773502691896257,
                                           0.5773502691896257};
constexpr std::array<double, 3> legendreX3{-0.7745966692414834, 0.0,
                                           0.7745966692414834};
constexpr std::array<double, 4> legendreX4{
    -0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
    0.8611363115940526};
// Legendre weighs
constexpr std::array<double, 1> legendreW1{2.0};
constexpr std::array<double, 2> legendreW2{1.0, 1.0};
constexpr std::array<double, 3> legendreW3{
    0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
constexpr std::array<double, 4> legendreW4{
    0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
    0.3478548451374538};

/* double legendreQuad2Do1(std::function<double(double x, double y)> f,
                        double xLow, double xHigh, double yLow, double yHigh) {

}

double legendreQuad2D(std::function<double(double x, double y)> f, double xLow,
                      double xHigh, double yLow, double yHigh, size_t order) {
  assert(order <= 4 && order > 0);
  switch (order) { case 1: }
} */

#endif
