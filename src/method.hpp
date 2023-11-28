#ifndef ALE_SOLVER_SRC_METHOD_HPP_
#define ALE_SOLVER_SRC_METHOD_HPP_

#include <barrier>
#include <future>

#include "grid.hpp"
#include "problem.hpp"

enum class MethodType { eALE };

enum class Status { eExit, eStepStart, eCalcdt };

class Method {
 public:
  Method(MethodType type, double tmin, double tmax, double CFL);
  Method(Method const &rhs) = default;
  Method(Method &&rhs) = default;

  Method &operator=(Method const &rhs) = delete;
  Method &operator=(Method &&rhs) = delete;

  virtual ~Method() = default;

 public:
  virtual void calc(double dt) = 0;
  virtual double calcdt() const = 0;

  virtual void dumpGrid(std::shared_ptr<Problem> problem) const = 0;

 public:
  double t;
  double CFL;

  const double tmin;
  const double tmax;

  const MethodType type;
};

class MethodALE : public Method {
 public:
  MethodALE(std::shared_ptr<Problem> problem, int sizeX, int sizeY, double epsx,
            double epsu, double CFL, int nThreads);

  virtual ~MethodALE();

 public:
  virtual void calc(double dt) override;
  virtual double calcdt() const override;

  virtual void dumpGrid(std::shared_ptr<Problem> problem) const override;

 private:
  void initializeMethod();
  void initializeThreads();

  void threadFunction(const int imin, const int imax, const int jmin,
                      const int jmax, const int threadNum);
  void parallelLagrangianPhase(const int imin, const int imax, const int jmin,
                               const int jmax, const int threadNum);
  void threadCalcdt(const int imin, const int imax, const int jmin,
                    const int jmax, const int threadNum) const;

  void resolveBoundaries();
  void resolveLeftBoundary();
  void resolveTopBoundary();
  void resolveRightBoundary();
  void resolveBottomBoundary();

  void checkPressureConvergence();
  void checkCoordsConvergence();

  void updateStatus(Status newStatus) const;

 private:
  std::shared_ptr<GridALE> grid;

  std::vector<std::vector<double>> x_;
  std::vector<std::vector<double>> y_;
  std::vector<std::vector<double>> xNext;
  std::vector<std::vector<double>> yNext;
  std::vector<std::vector<double>> uNext;
  std::vector<std::vector<double>> vNext;
  std::vector<std::vector<double>> rhoNext;
  std::vector<std::vector<double>> p_;
  std::vector<std::vector<double>> pNext;

  mutable std::vector<double> dts;
  double dt;

  const BoundaryType leftBoundary;
  const BoundaryType topBoundary;
  const BoundaryType rightBoundary;
  const BoundaryType bottomBoundary;

  mutable std::promise<Status> statusPromise;
  mutable std::shared_future<Status> statusFutureGlobal;

  std::vector<double> coordsConvergenceMax;
  std::vector<double> pressureConvergenceMax;

  bool bPressureConverged = false;
  bool bCoordsConverged = false;

  mutable std::barrier<> statusSync;
  std::barrier<> stepSync;
  mutable std::barrier<> calcdtSync;
  std::barrier<std::function<void()>> lagrangianBoundarySync;
  std::barrier<> lagrangianCoordsSync;
  std::barrier<std::function<void()>> lagrangianPressureConvergenceSync;
  std::barrier<std::function<void()>> lagrangianCoordsConvergenceSync;
  std::barrier<> lagrangianUpdateNodesSync;

  double epsx;
  double epsu;

  const double eps1 = 1e-4;
  const double eps2 = 1e-6;

  std::vector<std::future<void>> threads;
  const int nThreads;

  Logger logger;
};

#endif
