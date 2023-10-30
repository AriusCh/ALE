#include "method.hpp"

Method::Method(MethodType type) : type(type) {}
MethodType Method::getType() const { return type; }

MethodALE::MethodALE(std::unique_ptr<Problem> problem)
    : Method(MethodType::eALE), grids(problem->createALEGrids(100, 100)) {
  auto leftBoundaryTypes = problem->getLeftBoundaryTypes();
  auto topBoundaryTypes = problem->getTopBoundaryTypes();
  auto rightBoundaryTypes = problem->getRightBoundaryTypes();
  auto bottomBoundaryTypes = problem->getBottomBoundaryTypes();
  for (int i = 0; i < grids.size(); i++) {
    // left boundary
    boundaries.push_back(grids[i]->createExternalBoundary(
        leftBoundaryTypes[i], ExternalBoundarySide::eLeft));
    boundaries.push_back(grids[i]->createExternalBoundary(
        topBoundaryTypes[i], ExternalBoundarySide::eTop));
    boundaries.push_back(grids[i]->createExternalBoundary(
        rightBoundaryTypes[i], ExternalBoundarySide::eRight));
    boundaries.push_back(grids[i]->createExternalBoundary(
        bottomBoundaryTypes[i], ExternalBoundarySide::eBottom));
  }
}
void MethodALE::calc(double dt) {
  double epsx = 0.5, epsu = 0.5;

  std::vector<std::vector<std::vector<double>>> x_;
  std::vector<std::vector<std::vector<double>>> y_;
  std::vector<std::vector<std::vector<double>>> p_;

  for (auto &grid : grids) {
    x_.push_back(grid->x);
    y_.push_back(grid->y);
    p_.push_back(grid->p);
  }
  for (int num = 0; num < grids.size(); num++) {
    auto &grid = grids[num];

    for (int i = 1; i < grid->u.size() - 1; i++) {
      for (int j = 1; j < grid->u[0].size() - 1; j++) {
        double Fx =
            0.5 *
            (grid->p[i][j] * (grid->y[i + 1][j] - grid->y[i][j + 1]) +
             grid->p[i - 1][j] * (grid->y[i][j + 1] - grid->y[i - 1][j]) +
             grid->p[i - 1][j - 1] * (grid->y[i - 1][j] - grid->y[i][j - 1]) +
             grid->p[i][j - 1] * (grid->y[i][j - 1] - grid->y[i + 1][j]));
        double Fy =
            0.5 *
            (grid->p[i][j] * (grid->x[i][j + 1] - grid->x[i + 1][j]) +
             grid->p[i - 1][j] * (grid->x[i - 1][j] - grid->x[i][j + 1]) +
             grid->p[i - 1][j - 1] * (grid->x[i][j - 1] - grid->x[i - 1][j]) +
             grid->p[i][j - 1] * (grid->x[i + 1][j] - grid->x[i][j - 1]));
        double FxNext =
            0.5 *
            (p_[num][i][j] * (y_[num][i + 1][j] - y_[num][i][j + 1]) +
             p_[num][i - 1][j] * (y_[num][i][j + 1] - y_[num][i - 1][j]) +
             p_[num][i - 1][j - 1] * (y_[num][i - 1][j] - y_[num][i][j - 1]) +
             p_[num][i][j - 1] * (y_[num][i][j - 1] - y_[num][i + 1][j]));
        double FyNext =
            0.5 *
            (p_[num][i][j] * (x_[num][i][j + 1] - x_[num][i + 1][j]) +
             p_[num][i - 1][j] * (x_[num][i - 1][j] - x_[num][i][j + 1]) +
             p_[num][i - 1][j - 1] * (x_[num][i][j - 1] - x_[num][i - 1][j]) +
             p_[num][i][j - 1] * (x_[num][i + 1][j] - x_[num][i][j - 1]));

        grid->uNext[i][j] =
            grid->u[i][j] +
            dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
        grid->vNext[i][j] =
            grid->v[i][j] +
            dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
      }
    }
  }

  resolveBoundaries();

  for (int num = 0; num < grids.size(); num++) {
    auto &grid = grids[num];
    std::vector<std::vector<double>> xNext(grid->x);
    std::vector<std::vector<double>> yNext(grid->y);

    std::vector<std::vector<double>> volumeNext(
        grid->sizeX, std::vector<double>(grid->sizeY, 0.));

    for (int i = 0; i < grid->sizeX; i++) {
      for (int j = 0; j < grid->sizeY; j++) {
        xNext[i][j] = grid->x[i][j] + dt * (epsx * grid->u[i][j] +
                                            (1 - epsx) * grid->uNext[i][j]);
        yNext[i][j] = grid->y[i][j] + dt * (epsx * grid->v[i][j] +
                                            (1 - epsx) * grid->vNext[i][j]);

        volumeNext[i][j] = grid->getV(i, j, xNext, yNext);
        double rhoNext = grid->rho[i][j] * grid->getV(i, j, grid->x, grid->y) /
                         volumeNext[i][j];
        double e = grid->eos->gete(grid->rho[i][j], grid->p[i][j]);
        double eNext =
            e + grid->p[i][j] / grid->rho[i][j] *
                    (1 - volumeNext[i][j] / grid->getV(i, j, grid->x, grid->y));
        double pNext = grid->eos->getp(rhoNext, eNext);
        grid->rho[i][j] = rhoNext;
      }
    }
    grid->x = xNext;
    grid->y = xNext;
    grid->u = grid->uNext;
    grid->v = grid->vNext;
  }
}
double MethodALE::calcdt() const { return 0.002; }
void MethodALE::resolveBoundaries() {
  for (auto &boundary : boundaries) {
    boundary->resolveLagrangianPhase();
  }
}
