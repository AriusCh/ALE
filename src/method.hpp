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
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadKinematicCellMass(
      size_t celli, size_t cellj);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadThermoCellMass(
      size_t celli, size_t cellj);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadCellVectorMass();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  quadCellVectorLaplacian();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadCellMv(
      size_t celli, size_t cellj,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadCellM(
      size_t celli, size_t cellj,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadCellKv(
      size_t celli, size_t cellj,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadCellK(
      size_t celli, size_t cellj,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void quadVerticalFacesK(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void quadHorizontalFacesK(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
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
  size_t getAdvectionIndexFromCell(const size_t celli, const size_t cellj,
                                   const size_t k) const;
  size_t getQuadIndexFromCell(const size_t celli, const size_t cellj,
                              const size_t k) const;
  double getMinRhoFromCell(size_t celli, size_t cellj);
  double getMinRhoFromCellAndNeighbours(size_t celli, size_t cellj);
  double getMaxRhoFromCell(size_t celli, size_t cellj);
  double getMaxRhoFromCellAndNeighbours(size_t celli, size_t cellj);
  double getMinEFromCell(size_t celli, size_t cellj);
  double getMinEFromCellAndNeighbours(size_t celli, size_t cellj);
  double getMaxEFromCell(size_t celli, size_t cellj);
  double getMaxEFromCellAndNeighbours(size_t celli, size_t cellj);
  double getMinRhoFromRemapCell(size_t celli, size_t cellj);
  double getMinRhoFromRemapCellAndNeighbours(size_t celli, size_t cellj);
  double getMaxRhoFromRemapCell(size_t celli, size_t cellj);
  double getMaxRhoFromRemapCellAndNeighbours(size_t celli, size_t cellj);
  double getMinEFromRemapCell(size_t celli, size_t cellj);
  double getMinEFromRemapCellAndNeighbours(size_t celli, size_t cellj);
  double getMaxEFromRemapCell(size_t celli, size_t cellj);
  double getMaxEFromRemapCellAndNeighbours(size_t celli, size_t cellj);

  void initInitializers();
  void initBasisValues();
  void initQuadBasisValues();
  void initOutputBasisValues();
  void initRhoValues();
  void initKinematicVectors();
  void initThermodynamicVector();
  void initKinematicMassMatrix();
  void initThermodynamicInverseMassMatrix();
  void initForceMatrices();
  void initVectorMatrices();
  void initRemapMatrices();
  void initSolvers();

  void calcKinematicMassMatrix();
  void calcThermoMassMatrix();
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
  void calcRemapMatrices(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void calcMv(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
              const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void calcM(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
             const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void calcKv(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
              const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void calcK(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
             const Eigen::Matrix<double, Eigen::Dynamic, 1> &y);
  void calcAuxMat();

  void optimizeMesh();

  void transitionToRemap();
  void preTransitionRhoToRemap();
  void transitionRhoToRemap();
  void preRhoToRemap(size_t celli, size_t cellj);
  void postRhoToRemap(size_t celli, size_t cellj);
  void transitionEToRemap();
  void preEToRemap(size_t celli, size_t cellj);
  void postEToRemap(size_t celli, size_t cellj);

  void transitionToLagrange();
  void transitionRhoToLagrange();
  void transitionEToLagrange();
  void preEToLagrange(size_t celli, size_t cellj);
  void postEToLagrange(size_t celli, size_t cellj);

  void calcTau(double hmin, double soundSpeed, double rhoLocal,
               double maxViscosityCoeff);
  void calcRemapTau();

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

  void resolveBoundaryVector();
  void resolveLeftBoundaryVector();
  void resolveTopBoundaryVector();
  void resolveRightBoundaryVector();
  void resolveBottomBoundaryVector();

  void RK2step();
  void remap();
  void remapStep();

  void remapStepRhoPrepare();
  void remapStepRhoMinMaxCell(const size_t celli, const size_t cellj);
  void remapStepRhoAvgCell(const size_t celli, const size_t cellj);
  void remapStepRhoFlux();
  void remapStepRhoFluxLimiting();
  void remapStepRhoFinal();

  void remapStepRhoEPrepare();
  void remapStepRhoEMinMaxCell(const size_t celli, const size_t cellj);
  void remapStepRhoEAvgCell(const size_t celli, const size_t cellj);
  void remapStepRhoEFlux();
  void remapStepRhoEFluxLimiting();
  void remapStepRhoEFinal();

  // AUXILIARY FUNCTIONS
  Eigen::Matrix<double, 2, 2> getCellJacobian(
      const size_t celli, const size_t cellj, const size_t quadi,
      const size_t quadj, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) const;
  void FCTProjection();

 private:
  const size_t xSize;  // Number of cells in Ox direction
  const size_t ySize;  // Number of cells in Oy direction
  const size_t order;  // FEM order
  const size_t quadOrder;
  const size_t Nk;  // Number kinematic space points
  const size_t Nt;  // Number of thermodynamic space points
  const size_t Na;  // Number of advection space points
  const size_t remapFrequency = 1;
  const double q1 = 0.5;
  const double q2 = 2.0;
  const double alpha = 0.5;
  const double alphamu = 2.5;
  const double beta1 = 0.85;
  const double beta2 = 1.02;
  const double gamma = 0.8;
  const double eps = 1e-4;
  const double remapCFL = 0.9;
  const double remapTMax = 1.0;
  const double remapTMin = 0.0;
  double l0;

  Eigen::Matrix<double, Eigen::Dynamic, 1> x;
  Eigen::Matrix<double, Eigen::Dynamic, 1> y;
  Eigen::Matrix<double, Eigen::Dynamic, 1> u;
  Eigen::Matrix<double, Eigen::Dynamic, 1> v;
  Eigen::Matrix<double, Eigen::Dynamic, 1> e;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoQuad;

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

  Eigen::Matrix<double, Eigen::Dynamic, 1> xOptimal;
  Eigen::Matrix<double, Eigen::Dynamic, 1> yOptimal;
  Eigen::Matrix<double, Eigen::Dynamic, 1> uMesh;
  Eigen::Matrix<double, Eigen::Dynamic, 1> vMesh;

  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemap;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapHigh;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapLow;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapAvg;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapMin;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapMax;

  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoERemap;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoERemapHigh;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoERemapLow;
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoERemapAvg;
  Eigen::Matrix<double, Eigen::Dynamic, 1> eRemapMin;
  Eigen::Matrix<double, Eigen::Dynamic, 1> eRemapMax;

  Eigen::SparseMatrix<double> Mkx;
  Eigen::SparseMatrix<double> Mky;
  Eigen::SparseMatrix<double> Mt_inv;
  Eigen::SparseMatrix<double> Fx;
  Eigen::SparseMatrix<double> Fy;

  Eigen::SparseMatrix<double> vectorLaplacianX;
  Eigen::SparseMatrix<double> vectorLaplacianY;
  Eigen::SparseMatrix<double> vectorMassX;
  Eigen::SparseMatrix<double> vectorMassY;
  Eigen::SparseMatrix<double> Mv;
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Kv;
  Eigen::SparseMatrix<double> K;
  Eigen::SparseMatrix<double> Mlumped;
  Eigen::SparseMatrix<double> MlumpedInv;
  Eigen::SparseMatrix<double> MInv;
  Eigen::SparseMatrix<double> L;
  Eigen::SparseMatrix<double> Kupwinded;
  Eigen::SparseMatrix<double> D;
  Eigen::SparseMatrix<double> antidiffusiveFlux;
  Eigen::SparseMatrix<double> fluxLimitingFactors;

  // FCT PROJECTION
  // INPUT
  Eigen::Matrix<double, Eigen::Dynamic, 1> wQuad;
  Eigen::Matrix<double, Eigen::Dynamic, 1> detJacQuad;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> basisQuadValues;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xQuad;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xProj;
  Eigen::Matrix<double, Eigen::Dynamic, 1> yQuad;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ymin;
  Eigen::Matrix<double, Eigen::Dynamic, 1> ymax;
  // INTERNAL
  Eigen::Matrix<double, Eigen::Dynamic, 1> a;
  Eigen::Matrix<double, Eigen::Dynamic, 1> b;
  Eigen::Matrix<double, Eigen::Dynamic, 1> z;
  Eigen::Matrix<double, Eigen::Dynamic, 1> FCTCoeffPlus;
  Eigen::Matrix<double, Eigen::Dynamic, 1> FCTCoeffMinus;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> localMassMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> localMassMatrix_inv;
  Eigen::SparseMatrix<double> localMassMatrixLumped;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> localFluxMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> FCTCoeffMatrix;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xyHigh;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xyLow;
  // OUTPUT
  Eigen::Matrix<double, Eigen::Dynamic, 1> xyOut;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      kinematicSolverx;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      kinematicSolvery;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      MvSolver;

  size_t iteration = 0;
  double tau;
  double remapTau;
  double remapDt;
  double remapT;

  std::vector<double> kinematicBasis1DQuadValues;
  std::vector<double> kinematicBasis1DdxQuadValues;
  std::vector<double> thermoBasis1DQuadValues;
  std::vector<double> advectionBasis1DQuadValues;
  std::vector<double> advectionBasis1DdxQuadValues;

  std::vector<double> output1DKinematicValues;
  std::vector<double> output1DdxKinematicValues;
  std::vector<double> output1DThermoValues;
  std::vector<double> output1DQuadValues;

  std::function<std::shared_ptr<EOS>(double x, double y)>
      eosInitializer;  // Function that returns EOSes to use in the grid
};

#endif
