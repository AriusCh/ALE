#include "method.hpp"

#include <cassert>
#include <cmath>

Method::Method(MethodType type, double tmin, double tmax, double CFL)
    : t(tmin), CFL(CFL), tmin(tmin), tmax(tmax), type(type) {}

MethodALE::MethodALE(std::shared_ptr<Problem> problem, int sizeX, int sizeY,
                     double CFL, int nThreads)
    : Method(MethodType::eALE, problem->tmin, problem->tmax, CFL),
      grid(problem->createALEGrid(sizeX, sizeY)),
      leftBoundary(problem->leftBoundaryType),
      topBoundary(problem->topBoundaryType),
      rightBoundary(problem->rightBoundaryType),
      bottomBoundary(problem->bottomBoundaryType),
      lagrangianBoundarySync(nThreads, [this]() { resolveBoundaries(); }),
      lagrangianCoordsSync(nThreads),
      lagrangianPressureConvergenceSync(
          nThreads, [this]() { checkPressureConvergence(); }),
      lagrangianCoordsConvergenceSync(nThreads,
                                      [this]() { checkCoordsConvergence(); }),
      nThreads(nThreads) {
  initializeMethod();
  initializeThreads();
}
MethodALE::~MethodALE() {
  updateStatus(Status::eExit);
  for (auto &future : threads) {
    future.get();
  }
}
void MethodALE::calc(double dt) { this->dt = dt; }
double MethodALE::calcdt() const { return 0.002; }
void MethodALE::initializeMethod() {}
void MethodALE::initializeThreads() {
  int ncols = 1, nrows = 1;

  if (std::pow(std::floor(std::sqrt(nThreads)), 2) == nThreads) {
    ncols = std::sqrt(nThreads);
    nrows = std::sqrt(nThreads);
  } else {
    ncols = nThreads;
  }
  assert(ncols * nrows == nThreads);

  double di = static_cast<double>(grid->sizeX + 1) / static_cast<double>(ncols);
  double dj = static_cast<double>(grid->sizeY + 1) / static_cast<double>(nrows);
  threads.reserve(nThreads);
  for (int i = 0; i < nThreads; i++) {
    int imin = di * (i / ncols);
    int imax = i / ncols != ncols - 1 ? di * (i / ncols + 1) : grid->sizeX + 1;
    int jmin = dj * (i % nrows);
    int jmax = i % nrows == nrows - 1 ? dj * (i % nrows + 1) : grid->sizeY + 1;
    threads.push_back(std::async(std::launch::async, &MethodALE::threadFunction,
                                 this, imin, imax, jmin, jmax, i));
  }
}
void MethodALE::threadFunction(const int imin, const int imax, const int jmin,
                               const int jmax, const int threadNum) {
  while (true) {
    std::shared_future<Status> statusFuture(statusFutureGlobal);
    switch (statusFuture.get()) {
      case Status::eExit:
        return;
      case Status::eStepStart:
        parallelLagrangianPhase(imin, imax, jmin, jmax, threadNum);
        break;
    }
  }
}
void MethodALE::parallelLagrangianPhase(const int imin, const int imax,
                                        const int jmin, const int jmax,
                                        const int threadNum) {
  // assign Next values to current ones
  for (int i = imin; i < imax; i++) {
    for (int j = jmin; j < jmax; j++) {
      xNext[i][j] = grid->x[i][j];
      yNext[i][j] = grid->y[i][j];
      uNext[i][j] = grid->u[i][j];
      vNext[i][j] = grid->v[i][j];
      pNext[i][j] = grid->p[i][j];
    }
  }
  // Helper limits
  const double iminL = std::max(1, imin);
  const double imaxL = std::min(grid->sizeX, imax);
  const double jminL = std::max(1, jmax);
  const double jmaxL = std::min(grid->sizeY, jmax);
  do {
    // Step 1
    // Populate x_ and y_ with most recent guess
    for (int i = imin; i < imax; i++) {
      for (int j = jmin; j < jmax; j++) {
        x_[i][j] = xNext[i][j];
        y_[i][j] = yNext[i][j];
      }
    }

    do {
      // Step 2
      // Populate p_ with most recent guess
      for (int i = imin; i < imax; i++) {
        for (int j = jmin; j < jmax; j++) {
          p_[i][j] = pNext[i][j];
        }
      }

      // Step 3
      // Calculate velocity components uNext, vNext
      for (int i = iminL; i < imaxL; i++) {
        for (int j = jminL; j < jmaxL; j++) {
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
              0.5 * (p_[i][j] * (y_[i + 1][j] - y_[i][j + 1]) +
                     p_[i - 1][j] * (y_[i][j + 1] - y_[i - 1][j]) +
                     p_[i - 1][j - 1] * (y_[i - 1][j] - y_[i][j - 1]) +
                     p_[i][j - 1] * (y_[i][j - 1] - y_[i + 1][j]));
          double FyNext =
              0.5 * (p_[i][j] * (x_[i][j + 1] - x_[i + 1][j]) +
                     p_[i - 1][j] * (x_[i - 1][j] - x_[i][j + 1]) +
                     p_[i - 1][j - 1] * (x_[i][j - 1] - x_[i - 1][j]) +
                     p_[i][j - 1] * (x_[i + 1][j] - x_[i][j - 1]));

          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        }
      }
      // Wait for main thread to finish setting boundary conditions
      lagrangianBoundarySync.arrive_and_wait();

      // Step 4
      // Calculate coordinates of vertex (xNext, yNext)
      for (int i = imin; i < imax; i++) {
        for (int j = jmin; j < jmax; j++) {
          xNext[i][j] = grid->x[i][j] + dt * (epsx * grid->u[i][j] +
                                              (1.0 - epsx) * uNext[i][j]);
          yNext[i][j] = grid->y[i][j] + dt * (epsx * grid->v[i][j] +
                                              (1.0 - epsx) * vNext[i][j]);
        }
      }

      // Wait for all threads
      lagrangianCoordsSync.arrive_and_wait();

      // Step 5
      // Update pressures pNext
      pressureConvergenceMax[threadNum] = 0.0;
      for (int i = imin; i < imaxL; i++) {
        for (int j = jmin; j < jmaxL; j++) {
          double V = grid->getV(i, j, grid->x, grid->y);
          double V_ = grid->getV(i, j, xNext, yNext);

          double e = grid->eos->gete(grid->rho[i][j], grid->p[i][j]);
          double rhoNext = grid->rho[i][j] * V / V_;
          double eNext = e + grid->p[i][j] / grid->rho[i][j] * (1 - V_ / V);
          pNext[i][j] = grid->eos->getp(rhoNext, eNext);
          double dp = std::abs((p_[i][j] - pNext[i][j]) / p_[i][j]);
          if (dp > pressureConvergenceMax[threadNum]) {
            pressureConvergenceMax[threadNum] = dp;
          }
        }
      }

      // Sync and check the convergence of pressures
      lagrangianPressureConvergenceSync.arrive_and_wait();
    } while (!bPressureConverged);

    // Calc max dx
    coordsConvergenceMax[threadNum] = 0.0;
    for (int i = imin; i < imax; i++) {
      for (int j = jmin; j < jmax; j++) {
        double dx = std::abs((x_[i][j] - xNext[i][j]) / x_[i][j]);
        double dy = std::abs((y_[i][j] - yNext[i][j]) / y_[i][j]);
        if (std::max(dx, dy) > coordsConvergenceMax[threadNum]) {
          coordsConvergenceMax[threadNum] = std::max(dx, dy);
        }
      }
    }

    // Sync and check the convergence of coordinates
    lagrangianCoordsConvergenceSync.arrive_and_wait();
  } while (!bCoordsConverged);
}
void MethodALE::resolveBoundaries() {
  resolveLeftBoundary();
  resolveTopBoundary();
  resolveRightBoundary();
  resolveBottomBoundary();
}
void MethodALE::resolveLeftBoundary() {
  int i = 0;
  // resolve boundary then bottom left node
  switch (leftBoundary) {
    case BoundaryType::eExternalFree:
      for (int j = 1; j < grid->sizeY; j++) {
        double Fx =
            0.5 * (grid->p[i][j] * (grid->y[i + 1][j] - grid->y[i][j + 1]) +
                   grid->p[i][j - 1] * (grid->y[i][j - 1] - grid->y[i + 1][j]));
        double Fy =
            0.5 * (grid->p[i][j] * (grid->x[i][j + 1] - grid->x[i + 1][j]) +
                   grid->p[i][j - 1] * (grid->x[i + 1][j] - grid->x[i][j - 1]));
        double FxNext = 0.5 * (p_[i][j] * (y_[i + 1][j] - y_[i][j + 1]) +
                               p_[i][j - 1] * (y_[i][j - 1] - y_[i + 1][j]));
        double FyNext = 0.5 * (p_[i][j] * (x_[i][j + 1] - x_[i + 1][j]) +
                               p_[i][j - 1] * (x_[i + 1][j] - x_[i][j - 1]));

        uNext[i][j] = grid->u[i][j] +
                      dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
        vNext[i][j] = grid->v[i][j] +
                      dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
      }
      switch (bottomBoundary) {
        case BoundaryType::eExternalFree: {
          int j = 0;
          double Fx =
              0.5 * (grid->p[i][j] * (grid->y[i + 1][j] - grid->y[i][j + 1]));
          double Fy =
              0.5 * (grid->p[i][j] * (grid->x[i][j + 1] - grid->x[i + 1][j]));
          double FxNext = 0.5 * (p_[i][j] * (y_[i + 1][j] - y_[i][j + 1]));
          double FyNext = 0.5 * (p_[i][j] * (x_[i][j + 1] - x_[i + 1][j]));

          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        } break;
        case BoundaryType::eExternalTransparent: {
          int j = 0;
          double Fx =
              0.5 * (grid->p[i][j] * (grid->y[i + 1][j] - grid->y[i][j + 1]));
          double FxNext = 0.5 * (p_[i][j] * (y_[i + 1][j] - y_[i][j + 1]));
          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = vNext[i][j + 1];
        } break;
        case BoundaryType::eExternalNoSlipWall: {
          int j = 0;
          uNext[i][j] = 0.;
          vNext[i][j] = 0.;
        }
      }
      break;

    case BoundaryType::eExternalTransparent:
      for (int j = 1; j < grid->sizeY; j++) {
        uNext[i][j] = uNext[i + 1][j];
        vNext[i][j] = vNext[i + 1][j];
      }
      switch (bottomBoundary) {
        case BoundaryType::eExternalFree: {
          int j = 0;
          double Fy =
              0.5 * (grid->p[i][j] * (grid->x[i][j + 1] - grid->x[i + 1][j]));
          double FyNext = 0.5 * (p_[i][j] * (x_[i][j + 1] - x_[i + 1][j]));
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
          double Fx =
              0.5 *
              (grid->p[i + 1][j] * (grid->y[i + 2][j] - grid->y[i + 1][j + 1]) +
               grid->p[i][j] * (grid->y[i + 1][j + 1] - grid->y[i][j]));
          double FxNext =
              0.5 * (p_[i + 1][j] * (y_[i + 2][j] - y_[i + 1][j + 1]) +
                     p_[i][j] * (y_[i + 1][j + 1] - y_[i][j]));
          uNext[i][j] =
              grid->u[i + 1][j] +
              dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i + 1][j];
        } break;
        case BoundaryType::eExternalTransparent: {
          int j = 0;
          uNext[i][j] = uNext[i + 1][j + 1];
          vNext[i][j] = vNext[i + 1][j + 1];
        } break;
        case BoundaryType::eExternalNoSlipWall: {
          int j = 0;
          uNext[i][j] = 0.;
          vNext[i][j] = 0.;
        }
      }
      break;

    case BoundaryType::eExternalNoSlipWall:
      for (int j = 0; j < grid->sizeY; j++) {
        uNext[i][j] = 0.;
        vNext[i][j] = 0.;
      }
      break;
  }
}
void MethodALE::resolveTopBoundary() {
  int j = grid->sizeY;
  // Resolve top boundary and top left node
  switch (topBoundary) {
    case BoundaryType::eExternalFree:
      for (int i = 1; i < grid->sizeX; i++) {
        double Fx =
            0.5 *
            (grid->p[i - 1][j - 1] * (grid->y[i - 1][j] - grid->y[i][j - 1]) +
             grid->p[i][j - 1] * (grid->y[i][j - 1] - grid->y[i + 1][j]));
        double Fy =
            0.5 *
            (grid->p[i - 1][j - 1] * (grid->x[i][j - 1] - grid->x[i - 1][j]) +
             grid->p[i][j - 1] * (grid->x[i + 1][j] - grid->x[i][j - 1]));
        double FxNext =
            0.5 * (p_[i - 1][j - 1] * (y_[i - 1][j] - y_[i][j - 1]) +
                   p_[i][j - 1] * (y_[i][j - 1] - y_[i + 1][j]));
        double FyNext =
            0.5 * (p_[i - 1][j - 1] * (x_[i][j - 1] - x_[i - 1][j]) +
                   p_[i][j - 1] * (x_[i + 1][j] - x_[i][j - 1]));

        uNext[i][j] = grid->u[i][j] +
                      dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
        vNext[i][j] = grid->v[i][j] +
                      dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
      }
      switch (leftBoundary) {
        case BoundaryType::eExternalFree: {
          int i = 0;
          double Fx = 0.5 * (grid->p[i][j - 1] *
                             (grid->y[i][j - 1] - grid->y[i + 1][j]));
          double Fy = 0.5 * (grid->p[i][j - 1] *
                             (grid->x[i + 1][j] - grid->x[i][j - 1]));
          double FxNext = 0.5 * (p_[i][j - 1] * (y_[i][j - 1] - y_[i + 1][j]));
          double FyNext = 0.5 * (p_[i][j - 1] * (x_[i + 1][j] - x_[i][j - 1]));
          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        } break;
        case BoundaryType::eExternalTransparent: {
          int i = 0;
          double Fy = 0.5 * (grid->p[i][j - 1] *
                             (grid->x[i + 1][j] - grid->x[i][j - 1]));
          double FyNext = 0.5 * (p_[i][j - 1] * (x_[i + 1][j] - x_[i][j - 1]));
          uNext[i][j] = uNext[i + 1][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        } break;
        case BoundaryType::eExternalNoSlipWall: {
          int i = 0;
          uNext[i][j] = 0.;
          vNext[i][j] = 0.;
        }
      }
      break;
    case BoundaryType::eExternalTransparent:
      for (int i = 0; i < grid->sizeX; i++) {
        uNext[i][j] = uNext[i][j - 1];
        vNext[i][j] = vNext[i][j - 1];
      }
      break;
    case BoundaryType::eExternalNoSlipWall:
      for (int i = 0; i < grid->sizeX; i++) {
        uNext[i][j] = 0.;
        vNext[i][j] = 0.;
      }
      break;
  }
}
void MethodALE::resolveRightBoundary() {
  int i = grid->sizeX;
  switch (rightBoundary) {
    case BoundaryType::eExternalFree:
      for (int j = 1; j < grid->sizeY; j++) {
        double Fx =
            0.5 *
            (grid->p[i - 1][j] * (grid->y[i][j + 1] - grid->y[i - 1][j]) +
             grid->p[i - 1][j - 1] * (grid->y[i - 1][j] - grid->y[i][j - 1]));
        double Fy =
            0.5 *
            (grid->p[i - 1][j] * (grid->x[i - 1][j] - grid->x[i][j + 1]) +
             grid->p[i - 1][j - 1] * (grid->x[i][j - 1] - grid->x[i - 1][j]));
        double FxNext =
            0.5 * (p_[i - 1][j] * (y_[i][j + 1] - y_[i - 1][j]) +
                   p_[i - 1][j - 1] * (y_[i - 1][j] - y_[i][j - 1]));
        double FyNext =
            0.5 * (p_[i - 1][j] * (x_[i - 1][j] - x_[i][j + 1]) +
                   p_[i - 1][j - 1] * (x_[i][j - 1] - x_[i - 1][j]));

        uNext[i][j] = grid->u[i][j] +
                      dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
        vNext[i][j] = grid->v[i][j] +
                      dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
      }
      switch (topBoundary) {
        case BoundaryType::eExternalFree: {
          int j = grid->sizeY;
          double Fx = 0.5 * (grid->p[i - 1][j - 1] *
                             (grid->y[i - 1][j] - grid->y[i][j - 1]));
          double Fy = 0.5 * (grid->p[i - 1][j - 1] *
                             (grid->x[i][j - 1] - grid->x[i - 1][j]));
          double FxNext =
              0.5 * (p_[i - 1][j - 1] * (y_[i - 1][j] - y_[i][j - 1]));
          double FyNext =
              0.5 * (p_[i - 1][j - 1] * (x_[i][j - 1] - x_[i - 1][j]));

          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        } break;
        case BoundaryType::eExternalTransparent: {
          int j = grid->sizeY;
          uNext[i][j] = uNext[i][j - 1];
          vNext[i][j] = vNext[i][j - 1];
        } break;
        case BoundaryType::eExternalNoSlipWall: {
          int j = grid->sizeY;
          uNext[i][j] = 0.;
          vNext[i][j] = 0.;
        } break;
      }
      break;
    case BoundaryType::eExternalTransparent:
      for (int j = 1; j < grid->sizeY + 1; j++) {
        uNext[i][j] = uNext[i - 1][j];
        vNext[i][j] = vNext[i - 1][j];
      }
      break;
    case BoundaryType::eExternalNoSlipWall:
      for (int j = 1; j < grid->sizeY + 1; j++) {
        uNext[i][j] = 0.0;
        vNext[i][j] = 0.0;
      }
      break;
  }
}
void MethodALE::resolveBottomBoundary() {
  int j = 0;
  switch (bottomBoundary) {
    case BoundaryType::eExternalFree:
      for (int i = 1; i < grid->sizeX; i++) {
        double Fx =
            0.5 * (grid->p[i][j] * (grid->y[i + 1][j] - grid->y[i][j + 1]) +
                   grid->p[i - 1][j] * (grid->y[i][j + 1] - grid->y[i - 1][j]));
        double Fy =
            0.5 * (grid->p[i][j] * (grid->x[i][j + 1] - grid->x[i + 1][j]) +
                   grid->p[i - 1][j] * (grid->x[i - 1][j] - grid->x[i][j + 1]));
        double FxNext = 0.5 * (p_[i][j] * (y_[i + 1][j] - y_[i][j + 1]) +
                               p_[i - 1][j] * (y_[i][j + 1] - y_[i - 1][j]));
        double FyNext = 0.5 * (p_[i][j] * (x_[i][j + 1] - x_[i + 1][j]) +
                               p_[i - 1][j] * (x_[i - 1][j] - x_[i][j + 1]));

        uNext[i][j] = grid->u[i][j] +
                      dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
        vNext[i][j] = grid->v[i][j] +
                      dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
      }
      switch (rightBoundary) {
        case BoundaryType::eExternalFree: {
          int i = grid->sizeX;
          double Fx = 0.5 * (grid->p[i - 1][j] *
                             (grid->y[i][j + 1] - grid->y[i - 1][j]));
          double Fy = 0.5 * (grid->p[i - 1][j] *
                             (grid->x[i - 1][j] - grid->x[i][j + 1]));
          double FxNext = 0.5 * (p_[i - 1][j] * (y_[i][j + 1] - y_[i - 1][j]));
          double FyNext = 0.5 * (p_[i - 1][j] * (x_[i - 1][j] - x_[i][j + 1]));

          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        } break;
        case BoundaryType::eExternalTransparent: {
          int i = grid->sizeX;
          uNext[i][j] = uNext[i - 1][j];
          vNext[i][j] = vNext[i - 1][j];
        } break;
        case BoundaryType::eExternalNoSlipWall: {
          int i = grid->sizeX;
          uNext[i][j] = 0.0;
          vNext[i][j] = 0.0;

        } break;
      }
      break;
    case BoundaryType::eExternalTransparent:
      for (int i = 1; i < grid->sizeX + 1; i++) {
        uNext[i][j] = uNext[i][j + 1];
        vNext[i][j] = vNext[i][j + 1];
      }
      break;
    case BoundaryType::eExternalNoSlipWall:
      for (int i = 1; i < grid->sizeX + 1; i++) {
        uNext[i][j] = 0.0;
        vNext[i][j] = 0.0;
      }
      break;
  }
}
void MethodALE::checkPressureConvergence() {
  if (std::any_of(pressureConvergenceMax.begin(), pressureConvergenceMax.end(),
                  [eps = this->eps1](double dp) { return dp > eps; })) {
    bPressureConverged = false;
    return;
  }
  bPressureConverged = true;
}
void MethodALE::checkCoordsConvergence() {
  if (std::any_of(coordsConvergenceMax.begin(), pressureConvergenceMax.end(),
                  [eps = this->eps2](double dx) { return dx > eps; })) {
    bCoordsConverged = false;
    return;
  }
  bCoordsConverged = true;
}
void MethodALE::updateStatus(Status newStatus) {
  std::promise<Status> newPromise;
  statusFutureGlobal = newPromise.get_future();
  statusPromise.set_value(newStatus);
  statusPromise = std::move(newPromise);
}
