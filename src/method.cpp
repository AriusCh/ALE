#include "method.hpp"

#include <cassert>
#include <cmath>
#include <format>

Method::Method(MethodType type, double tmin, double tmax, double CFL)
    : t(tmin), CFL(CFL), tmin(tmin), tmax(tmax), type(type) {}

MethodALE::MethodALE(std::shared_ptr<Problem> problem, int sizeX, int sizeY,
                     double epsx, double epsu, double CFL, int nThreads)
    : Method(MethodType::eALE, problem->tmin, problem->tmax, CFL),
      grid(problem->createALEGrid(sizeX, sizeY)),
      leftBoundary(problem->leftBoundaryType),
      topBoundary(problem->topBoundaryType),
      rightBoundary(problem->rightBoundaryType),
      bottomBoundary(problem->bottomBoundaryType),
      coordsConvergenceMax(nThreads),
      pressureConvergenceMax(nThreads),
      statusSync(nThreads + 1),
      stepSync(nThreads + 1),
      lagrangianBoundarySync(nThreads, [this]() { resolveBoundaries(); }),
      lagrangianCoordsSync(nThreads),
      lagrangianPressureConvergenceSync(
          nThreads, [this]() { checkPressureConvergence(); }),
      lagrangianCoordsConvergenceSync(nThreads,
                                      [this]() { checkCoordsConvergence(); }),
      lagrangianUpdateNodesSync(nThreads),
      epsx(epsx),
      epsu(epsu),
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
void MethodALE::calc(double dt) {
  this->dt = dt;
  updateStatus(Status::eStepStart);
  stepSync.arrive_and_wait();
}
double MethodALE::calcdt() const {
  double dt = tmax - tmin;
  for (int i = 0; i < grid->sizeX; i++) {
    for (int j = 0; j < grid->sizeY; j++) {
      double dx = grid->x[i + 1][j] - grid->x[i][j];
      double dy = grid->y[i][j + 1] - grid->y[i][j];
      double u = 0.25 * (grid->u[i + 1][j] + grid->u[i + 1][j + 1] +
                         grid->u[i][j + 1] + grid->u[i][j]);
      double v = 0.25 * (grid->v[i + 1][j] + grid->v[i + 1][j + 1] +
                         grid->v[i][j + 1] + grid->v[i][j]);
      double c = grid->eos->getc(grid->rho[i][j], grid->p[i][j]);
      double dt1 = dx / (c + std::abs(u));
      double dt2 = dy / (c + std::abs(v));
      if (std::min(dt1, dt2) < dt) {
        dt = std::min(dt1, dt2);
      }
    }
  }
  return CFL * dt;

  return 0.001;
}
void MethodALE::dumpGrid(std::shared_ptr<Problem> problem) const {
  problem->dumpGrid(grid, t);
}
void MethodALE::initializeMethod() {
  x_ = grid->x;
  y_ = grid->y;
  xNext = grid->x;
  yNext = grid->y;
  uNext = grid->u;
  vNext = grid->v;
  rhoNext = grid->rho;
  p_ = grid->p;
  pNext = grid->p;

  statusFutureGlobal = statusPromise.get_future();
}
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
    int imin = di * (i / nrows);
    int imax = i / nrows != ncols - 1 ? di * (i / nrows + 1) : grid->sizeX + 1;
    int jmin = dj * (i % nrows);
    int jmax = i % nrows != nrows - 1 ? dj * (i % nrows + 1) : grid->sizeY + 1;
    threads.push_back(std::async(std::launch::async, &MethodALE::threadFunction,
                                 this, imin, imax, jmin, jmax, i));
  }
}
void MethodALE::threadFunction(const int imin, const int imax, const int jmin,
                               const int jmax, const int threadNum) {
  while (true) {
    std::shared_future<Status> statusFuture(statusFutureGlobal);
    statusSync.arrive_and_wait();
    switch (statusFuture.get()) {
      case Status::eExit:
        return;
      case Status::eStepStart:
        parallelLagrangianPhase(imin, imax, jmin, jmax, threadNum);
        stepSync.arrive_and_wait();
        break;
    }
  }
}
void MethodALE::parallelLagrangianPhase(const int imin, const int imax,
                                        const int jmin, const int jmax,
                                        const int threadNum) {
  // Helper limits
  const double iminL = std::max(1, imin);
  const double imaxL = std::min(grid->sizeX, imax);
  const double jminL = std::max(1, jmin);
  const double jmaxL = std::min(grid->sizeY, jmax);
  // assign Next values to current ones
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = imin; i < imax; i++) {
    for (int j = jmin; j < jmax; j++) {
      xNext[i][j] = grid->x[i][j] + dt * grid->u[i][j];
      yNext[i][j] = grid->y[i][j] + dt * grid->v[i][j];
      uNext[i][j] = grid->u[i][j];
      vNext[i][j] = grid->v[i][j];
    }
  }
  for (int i = imin; i < imaxL; i++) {
    for (int j = jmin; j < jmaxL; j++) {
      rhoNext[i][j] = grid->rho[i][j];
      pNext[i][j] = grid->p[i][j];
    }
  }
  if (threadNum == 0) {
    auto end = std::chrono::high_resolution_clock::now();
    std::string message =
        std::format("THREAD_{} ASSIGN TIME: {:.6f}", threadNum,
                    std::chrono::duration<double>(end - start).count());
    logger.Log(message);
  }
  do {
    // Step 1
    start = std::chrono::high_resolution_clock::now();
    // Populate x_ and y_ with most recent guess
    for (int i = imin; i < imax; i++) {
      for (int j = jmin; j < jmax; j++) {
        x_[i][j] = xNext[i][j];
        y_[i][j] = yNext[i][j];
      }
    }
    if (threadNum == 0) {
      auto end = std::chrono::high_resolution_clock::now();
      std::string message =
          std::format("THREAD_{} STEP 1 TIME: {:.6f}", threadNum,
                      std::chrono::duration<double>(end - start).count());
      logger.Log(message);
    }

    do {
      // Step 2
      start = std::chrono::high_resolution_clock::now();
      // Populate p_ with most recent guess
      for (int i = imin; i < imaxL; i++) {
        for (int j = jmin; j < jmaxL; j++) {
          p_[i][j] = pNext[i][j];
        }
      }
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message =
            std::format("THREAD_{} STEP 2 TIME: {:.6f}", threadNum,
                        std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }

      // Step 3
      start = std::chrono::high_resolution_clock::now();
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
              0.5 * (pNext[i][j] * (yNext[i + 1][j] - yNext[i][j + 1]) +
                     pNext[i - 1][j] * (yNext[i][j + 1] - yNext[i - 1][j]) +
                     pNext[i - 1][j - 1] * (yNext[i - 1][j] - yNext[i][j - 1]) +
                     pNext[i][j - 1] * (yNext[i][j - 1] - yNext[i + 1][j]));
          double FyNext =
              0.5 * (pNext[i][j] * (xNext[i][j + 1] - xNext[i + 1][j]) +
                     pNext[i - 1][j] * (xNext[i - 1][j] - xNext[i][j + 1]) +
                     pNext[i - 1][j - 1] * (xNext[i][j - 1] - xNext[i - 1][j]) +
                     pNext[i][j - 1] * (xNext[i + 1][j] - xNext[i][j - 1]));

          uNext[i][j] = grid->u[i][j] +
                        dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
          vNext[i][j] = grid->v[i][j] +
                        dt * (epsu * Fy + (1 - epsu) * FyNext) / grid->m[i][j];
        }
      }
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message =
            std::format("THREAD_{} STEP 3 TIME: {:.6f}", threadNum,
                        std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }
      // Wait for all threads to finish and set boundary conditions
      start = std::chrono::high_resolution_clock::now();
      lagrangianBoundarySync.arrive_and_wait();
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message =
            std::format("THREAD_{} BOUNDARY SYNC TIME: {:.6f}", threadNum,
                        std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }

      // Step 4
      // Calculate coordinates of vertex (xNext, yNext)
      start = std::chrono::high_resolution_clock::now();
      for (int i = imin; i < imax; i++) {
        for (int j = jmin; j < jmax; j++) {
          xNext[i][j] = grid->x[i][j] + dt * (epsx * grid->u[i][j] +
                                              (1.0 - epsx) * uNext[i][j]);
          yNext[i][j] = grid->y[i][j] + dt * (epsx * grid->v[i][j] +
                                              (1.0 - epsx) * vNext[i][j]);
        }
      }
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message =
            std::format("THREAD_{} STEP 4 TIME: {:.6f}", threadNum,
                        std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }

      // Wait for all threads
      start = std::chrono::high_resolution_clock::now();
      lagrangianCoordsSync.arrive_and_wait();
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message =
            std::format("THREAD_{} COORDS SYNC TIME: {:.6f}", threadNum,
                        std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }

      // Step 5
      // Update pressures pNext
      start = std::chrono::high_resolution_clock::now();
      pressureConvergenceMax[threadNum] = 0.0;
      for (int i = imin; i < imaxL; i++) {
        for (int j = jmin; j < jmaxL; j++) {
          double V = grid->getV(i, j, grid->x, grid->y);
          double V_ = grid->getV(i, j, xNext, yNext);

          double e = grid->eos->gete(grid->rho[i][j], grid->p[i][j]);
          rhoNext[i][j] = grid->rho[i][j] * V / V_;
          double eNext = e + grid->p[i][j] / grid->rho[i][j] * (1 - V_ / V);
          pNext[i][j] = grid->eos->getp(rhoNext[i][j], eNext);
          // pNext[i][j] = grid->eos->gets(grid->rho[i][j], grid->p[i][j]) *
          //               std::pow(rhoNext, 1.4);
          double dp = std::abs(p_[i][j] - pNext[i][j]);
          if (dp > pressureConvergenceMax[threadNum]) {
            pressureConvergenceMax[threadNum] = dp;
          }
        }
      }
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message =
            std::format("THREAD_{} STEP 5 TIME: {:.6f}", threadNum,
                        std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }

      // Sync and check the convergence of pressures
      start = std::chrono::high_resolution_clock::now();
      lagrangianPressureConvergenceSync.arrive_and_wait();
      if (threadNum == 0) {
        auto end = std::chrono::high_resolution_clock::now();
        std::string message = std::format(
            "THREAD_{} PRESSURE CONVERGENCE SYNC TIME: {:.6f}", threadNum,
            std::chrono::duration<double>(end - start).count());
        logger.Log(message);
      }
    } while (!bPressureConverged);

    // Calc max dx
    start = std::chrono::high_resolution_clock::now();
    coordsConvergenceMax[threadNum] = 0.0;
    for (int i = imin; i < imax; i++) {
      for (int j = jmin; j < jmax; j++) {
        double dx = std::abs(x_[i][j] - xNext[i][j]);
        double dy = std::abs(y_[i][j] - yNext[i][j]);
        if (std::max(dx, dy) > coordsConvergenceMax[threadNum]) {
          coordsConvergenceMax[threadNum] = std::max(dx, dy);
        }
      }
    }
    if (threadNum == 0) {
      auto end = std::chrono::high_resolution_clock::now();
      std::string message =
          std::format("THREAD_{} COORDS CONVERGENCE TIME: {:.6f}", threadNum,
                      std::chrono::duration<double>(end - start).count());
      logger.Log(message);
    }

    // Sync and check the convergence of coordinates
    start = std::chrono::high_resolution_clock::now();
    lagrangianCoordsConvergenceSync.arrive_and_wait();
    if (threadNum == 0) {
      auto end = std::chrono::high_resolution_clock::now();
      std::string message = std::format(
          "THREAD_{} COORDS CONVERGENCE SYNC TIME: {:.6f}", threadNum,
          std::chrono::duration<double>(end - start).count());
      logger.Log(message);
    }
  } while (!bCoordsConverged);

  // Update densities, pressures and total energy
  start = std::chrono::high_resolution_clock::now();
  for (int i = imin; i < imaxL; i++) {
    for (int j = jmin; j < jmaxL; j++) {
      grid->rho[i][j] = rhoNext[i][j];

      grid->p[i][j] = pNext[i][j];
    }
  }
  // Update coordinates and velocities
  for (int i = imin; i < imax; i++) {
    for (int j = jmin; j < jmax; j++) {
      grid->x[i][j] = xNext[i][j];
      grid->y[i][j] = yNext[i][j];
      grid->u[i][j] = uNext[i][j];
      grid->v[i][j] = vNext[i][j];
    }
  }
  if (threadNum == 0) {
    auto end = std::chrono::high_resolution_clock::now();
    std::string message =
        std::format("THREAD_{} UPDATE TIME: {:.6f}", threadNum,
                    std::chrono::duration<double>(end - start).count());
    logger.Log(message);
  }
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
  double maxp = 0.0;
  for (int i = 0; i < grid->sizeX; i++) {
    for (int j = 0; j < grid->sizeY; j++) {
      if (std::abs(grid->p[i][j]) > maxp) {
        maxp = std::abs(grid->p[i][j]);
      }
    }
  }
  if (std::any_of(pressureConvergenceMax.begin(), pressureConvergenceMax.end(),
                  [eps = this->eps1 * maxp](double dp) { return dp > eps; })) {
    bPressureConverged = false;
    return;
  }
  bPressureConverged = true;
}
void MethodALE::checkCoordsConvergence() {
  double maxx = 0.0;
  for (int i = 0; i < grid->sizeX + 1; i++) {
    for (int j = 0; j < grid->sizeY + 1; j++) {
      if (std::abs(grid->x[i][j]) > maxx) {
        maxx = std::abs(grid->x[i][j]);
      }
      if (std::abs(grid->y[i][j]) > maxx) {
        maxx = std::abs(grid->y[i][j]);
      }
    }
  }
  if (std::any_of(coordsConvergenceMax.begin(), coordsConvergenceMax.end(),
                  [eps = this->eps2 * maxx](double dx) { return dx > eps; })) {
    bCoordsConverged = false;
    return;
  }
  bCoordsConverged = true;
}
void MethodALE::updateStatus(Status newStatus) {
  statusSync.arrive_and_wait();
  std::promise<Status> newPromise;
  statusFutureGlobal = newPromise.get_future();
  statusPromise.set_value(newStatus);
  statusPromise = std::move(newPromise);
}
