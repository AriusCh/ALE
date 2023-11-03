#ifndef ALE_SOLVER_SRC_METHOD_HPP_
#define ALE_SOLVER_SRC_METHOD_HPP_

#include <barrier>
#include <future>

#include "grid.hpp"
#include "problem.hpp"

enum class MethodType { eALE };

enum class Status { eExit, eStepStart };

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

 public:
  double t;
  double CFL;

  const double tmin;
  const double tmax;

  const MethodType type;
};

class MethodALE : public Method {
 public:
  MethodALE(std::shared_ptr<Problem> problem, int sizeX, int sizeY, double CFL,
            int nThreads = 1);

 public:
  virtual void calc(double dt);
  virtual double calcdt() const;

 private:
  void initializeMethod();
  void initializeThreads();

  void resolveBoundaries();

  void checkPressureConvergence();
  void checkCoordsConvergence();

  void updateStatus(Status newStatus);

 private:
  std::shared_ptr<GridALE> grid;

  std::vector<std::vector<double>> x_;
  std::vector<std::vector<double>> y_;
  std::vector<std::vector<double>> xNext;
  std::vector<std::vector<double>> yNext;
  std::vector<std::vector<double>> uNext;
  std::vector<std::vector<double>> vNext;
  std::vector<std::vector<double>> p_;
  std::vector<std::vector<double>> pNext;

  double dt;

  const BoundaryType leftBoundary;
  const BoundaryType topBoundary;
  const BoundaryType rightBoundary;
  const BoundaryType bottomBoundary;

  std::promise<Status> statusPromise;
  std::shared_future<Status> statusFuture;

  std::vector<double> coordsConvergenceMax;
  std::vector<double> pressureConvergenceMax;

  bool bPressureConverged = false;
  bool bCoordsConverged = false;

  std::barrier<std::function<void()>> lagrangianBoundarySync;
  std::barrier<> lagrangianCoordsSync;
  std::barrier<std::function<void()>> lagrangianPressureConvergenceSync;
  std::barrier<std::function<void()>> lagrangianCoordsConvergenceSync;

  std::function<void(const int imin, const int imax, const int jmin,
                     const int jmax, const int threadNum)>
      threadFunction;
  std::function<void(const int imin, const int imax, const int jmin,
                     const int jmax, const int threadNum)>
      lagrangianPhaseFunction;

  double epsx = 0.5;
  double epsu = 0.5;

  const double eps1 = 1e-4;
  const double eps2 = 1e-6;

  std::vector<std::future<void>> threads;
  const int nThreads;
};

#endif
