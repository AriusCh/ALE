#ifndef ALE_SOLVER_SRC_METHOD_HPP_
#define ALE_SOLVER_SRC_METHOD_HPP_

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>

#include "output_mgr.hpp"
#include "problem.hpp"

enum class Status { eExit, eStepStart, eCalcdt };

class Method {
 public:
  Method(const std::string &name_, const Problem &problem_)
      : name(name_),
        problem(problem_),
        dt(problem.tmax - problem.tmin),
        t(problem.tmin),
        writer(problem.name, name) {}
  Method(Method const &rhs) = default;
  Method(Method &&rhs) = default;

  Method &operator=(Method const &rhs) = delete;
  Method &operator=(Method &&rhs) = delete;

  virtual ~Method() = default;

 public:
  virtual void calc() = 0;
  virtual void calcdt() const = 0;

  virtual void dumpData() const = 0;
  virtual void dumpGrid() const = 0;

 public:
  std::string name;
  Problem problem;

  double dt;
  double t;

  Writer writer;
};

class FEMALEMethod : public Method {
 public:
  FEMALEMethod(const std::string &name, const Problem &problem, size_t xSize_,
               size_t ySize_, size_t order_);
  FEMALEMethod(FEMALEMethod const &rhs) = delete;
  FEMALEMethod(FEMALEMethod &&rhs) = delete;

  FEMALEMethod &operator=(FEMALEMethod const &rhs) = delete;
  FEMALEMethod &operator=(FEMALEMethod &&rhs) = delete;

  virtual ~FEMALEMethod() = default;

 public:
  virtual void calc();
  virtual void calcdt() const;

  virtual void dumpData() const;
  virtual void dumpGrid() const;

 private:
  double quadKinematicCellMass(size_t celli, size_t cellj, size_t basisi,
                               size_t basisj);
  double quadThermoCellMass(size_t celli, size_t cellj, size_t basisi,
                            size_t basisj);
  std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 2>
  quadForceCell(size_t celli, size_t cellj,
                const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
                const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
                const Eigen::Matrix<double, Eigen::Dynamic, 1> &e);

  size_t getKinematicIndexFromCell(const size_t celli, const size_t cellj,
                                   const size_t k) const;
  size_t getThermodynamicIndexFromCell(const size_t celli, const size_t cellj,
                                       const size_t k) const;

  void initInitializers();
  void initBasisValues();
  void initKinematicMassBasisValues();
  void initThermodynamicMassBasisValues();
  void initForceBasisValues();
  void initOutputBasisValues();
  void initKinematicVectors();
  void initThermodynamicVector();
  void initKinematicMassMatrix();
  void initThermodynamicInverseMassMatrix();
  void initForceMatrices();
  void initSolvers();

  void calcForceMatrices(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &e);
  Eigen::Matrix<double, 2, 2> calcStressTensor(
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &e,
      const Eigen::Matrix<double, 2, 2> &jacobian, double &soundSpeed,
      double &rhoLocal, double &maxViscosityCoeff, const size_t celli,
      const size_t cellj, const size_t i, const size_t j) const;
  Eigen::Matrix<double, 2, 2> calcArtificialViscosity(
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
      const Eigen::Matrix<double, 2, 2> &jacobian,
      const Eigen::Matrix<double, 2, 2> &jacobianInitial, double soundSpeed,
      double rhoLocal, double &maxViscosityCoeff, const size_t celli,
      const size_t cellj, const size_t i, const size_t j) const;
  double calcViscosityCoeff(const Eigen::Matrix<double, 2, 2> &jacobian,
                            const Eigen::Matrix<double, 2, 2> &jacobianInitial,
                            const Eigen::Matrix<double, 2, 2> &velocityGrad,
                            double eigenvalue,
                            const Eigen::Matrix<double, 2, 1> &eigenvector,
                            double soundSpeed, double rhoLocal) const;

  void calcTau(double hmin, double soundSpeed, double rhoLocal,
               double maxViscosityCoeff);

  void resolveBoundaryMass();
  void resolveLeftBoundaryMass();
  void resolveTopBoundaryMass();
  void resolveRightBoundaryMass();
  void resolveBottomBoundaryMass();
  void resolveBoundaryForce();
  void resolveLeftBoundaryForce();
  void resolveTopBoundaryForce();
  void resolveRightBoundaryForce();
  void resolveBottomBoundaryForce();

  void RK2step();

 private:
  const size_t xSize;  // Number of cells in Ox direction
  const size_t ySize;  // Number of cells in Oy direction
  const size_t order;  // FEM order
  const size_t kinematicMassQuadOrder;
  const size_t thermoMassQuadOrder;
  const size_t forceQuadOrder;
  const size_t Nk;  // Number kinematic space points
  const size_t Nt;  // Number of thermodynamic space points
  const double q1 = 0.5;
  const double q2 = 2.0;
  const double alpha = 0.5;
  const double alphamu = 2.5;
  const double beta1 = 0.85;
  const double beta2 = 1.02;
  const double gamma = 0.8;
  double l0;

  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  Eigen::Matrix<double, Eigen::Dynamic, 1> u;
  Eigen::Matrix<double, Eigen::Dynamic, 1> v;
  Eigen::Matrix<double, Eigen::Dynamic, 1> e;

  Eigen::Matrix<double, Eigen::Dynamic, 1> xInitial;
  Eigen::Matrix<double, Eigen::Dynamic, 1> yInitial;

  Eigen::Matrix<double, Eigen::Dynamic, 1> x05;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y05;
  Eigen::Matrix<double, Eigen::Dynamic, 1> u05;
  Eigen::Matrix<double, Eigen::Dynamic, 1> v05;
  Eigen::Matrix<double, Eigen::Dynamic, 1> e05;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Fu;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Fv;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Fe;

  Eigen::SparseMatrix<double> Mkx;
  Eigen::SparseMatrix<double> Mky;
  Eigen::SparseMatrix<double> Mt_inv;
  Eigen::SparseMatrix<double> Fx;
  Eigen::SparseMatrix<double> Fy;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      kinematicSolverx;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      kinematicSolvery;

  double tau;

  std::vector<double> kinematicMass1DValues;
  std::vector<double> kinematicMass1DdxValues;

  std::vector<double> thermoMass1DKinematicValues;
  std::vector<double> thermoMass1DdxKinematicValues;
  std::vector<double> thermoMass1DThermoValues;

  std::vector<double> force1DKinematicValues;
  std::vector<double> force1DdxKinematicValues;
  std::vector<double> force1DThermoValues;

  std::vector<double> output1DKinematicValues;
  std::vector<double> output1DdxKinematicValues;
  std::vector<double> output1DThermoValues;

  std::function<double(double x, double y)>
      rhoInitializer;  // Function that returns rho initial values
  std::function<std::shared_ptr<EOS>(double x, double y)>
      eosInitializer;  // Function that returns EOSes to use in the grid
};

#endif
