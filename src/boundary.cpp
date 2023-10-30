#include "boundary.hpp"
Boundary::Boundary(BoundaryType type_) : type(type_) {}

ExternalBoundary::ExternalBoundary(BoundaryType type_,
                                   ExternalBoundarySide side_,
                                   std::vector<std::vector<double>> &x,
                                   std::vector<std::vector<double>> &y,
                                   std::vector<std::vector<double>> &rho,
                                   std::vector<std::vector<double>> &p,
                                   std::vector<std::vector<double>> &m,
                                   std::vector<std::vector<double>> &uNext,
                                   std::vector<std::vector<double>> &vNext)
    : Boundary(type_),
      side(side_),
      x(x),
      y(y),
      rho(rho),
      p(p),
      m(m),
      uNext(uNext),
      vNext(vNext) {}
void ExternalBoundary::resolveLagrangianPhase() {
  switch (getType()) {
    case BoundaryType::eExternalTransparent:
      resolveLagrangianPhaseExternalTransparent(uNext, vNext);
      break;
  }
}
void ExternalBoundary::resolveLagrangianPhaseExternalTransparent(
    std::vector<std::vector<double>> &uNext,
    std::vector<std::vector<double>> &vNext) {
  switch (side) {
    case ExternalBoundarySide::eLeft:
      for (int j = 0; j < uNext[0].size(); j++) {
        const int i = 0;
        uNext[i][j] = uNext[i + 1][j];
        vNext[i][j] = vNext[i + 1][j];
      }
      break;
    case ExternalBoundarySide::eTop:
      for (int i = 0; i < uNext.size(); i++) {
        const int j = uNext[0].size() - 1;
        uNext[i][j] = uNext[i][j - 1];
        vNext[i][j] = vNext[i][j - 1];
      }
      break;
    case ExternalBoundarySide::eRight:
      for (int j = 0; j < uNext[0].size(); j++) {
        const int i = uNext.size() - 1;
        uNext[i][j] = uNext[i - 1][j];
        vNext[i][j] = vNext[i - 1][j];
      }
      break;
    case ExternalBoundarySide::eBottom:
      for (int i = 0; i < uNext.size(); i++) {
        const int j = 0;
        uNext[i][j] = uNext[i][j + 1];
        vNext[i][j] = vNext[i][j + 1];
      }
      break;
  }
}
