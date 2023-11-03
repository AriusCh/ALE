#include "method.hpp"

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
void MethodALE::calc(double dt) { this->dt = dt; }
double MethodALE::calcdt() const { return 0.002; }
void MethodALE::initializeMethod() {
  threadFunction = [&statusFutureGlobal = this->statusFuture,
                    lagrangianPhaseFunction = this->lagrangianPhaseFunction](
                       const int imin, const int imax, const int jmin,
                       const int jmax, const int threadNum) {
    while (true) {
      std::shared_future<Status> statusFuture(statusFutureGlobal);
      switch (statusFuture.get()) {
        case Status::eExit:
          return;
        case Status::eStepStart:
          lagrangianPhaseFunction(imin, imax, jmin, jmax, threadNum);
          break;
      }
    }
  };
  lagrangianPhaseFunction = [this](const int imin, const int imax,
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
                 grid->p[i - 1][j - 1] *
                     (grid->y[i - 1][j] - grid->y[i][j - 1]) +
                 grid->p[i][j - 1] * (grid->y[i][j - 1] - grid->y[i + 1][j]));
            double Fy =
                0.5 *
                (grid->p[i][j] * (grid->x[i][j + 1] - grid->x[i + 1][j]) +
                 grid->p[i - 1][j] * (grid->x[i - 1][j] - grid->x[i][j + 1]) +
                 grid->p[i - 1][j - 1] *
                     (grid->x[i][j - 1] - grid->x[i - 1][j]) +
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

            uNext[i][j] =
                grid->u[i][j] +
                dt * (epsu * Fx + (1 - epsu) * FxNext) / grid->m[i][j];
            vNext[i][j] =
                grid->v[i][j] +
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
  };
}
void MethodALE::initializeThreads() {
  threads.reserve(nThreads);
  for (int i = 0; i < nThreads; i++) {
    threads.push_back(std::async(threadFunction, 0, 0, 0, 0, i));
  }
}
void MethodALE::resolveBoundaries() {}
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
  statusFuture = newPromise.get_future();
  statusPromise.set_value(newStatus);
  statusPromise = std::move(newPromise);
}
