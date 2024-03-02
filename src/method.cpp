#include "method.hpp"

#include <algorithm>
#include <cassert>
#include <format>

#include "nodes.hpp"

FEMALEMethod::FEMALEMethod(const std::string &name, const Problem &problem_,
                           size_t xSize_, size_t ySize_, size_t order_)
    :  // Method(name, pr_->name, pr_->tmin, pr_->tmax, pr_->tOut, pr_->tMul),
      Method(name, problem_),
      xSize(xSize_),
      ySize(ySize_),
      order(order_),
      kinematicMassQuadOrder(order * 2),
      thermoMassQuadOrder(order * 2),
      forceQuadOrder(order * 2),
      Nk((order * xSize + 1) * (order * ySize + 1)),
      Nt((order * xSize) * (order * ySize)),
      // xmin(problem_->xmin),
      // xmax(problem_->xmax),
      // ymin(problem_->ymin),
      // ymax(problem_->ymax),
      // leftBoundaryType(problem_->leftBoundaryType),
      // topBoundaryType(problem_->topBoundaryType),
      // rightBoundaryType(problem_->rightBoundaryType),
      // bottomBoundaryType(problem_->bottomBoundaryType),
      // dimension(problem_->dimension),
      x(Nk),
      y(Nk),
      u(Nk),
      v(Nk),
      e(Nt),
      xInitial(Nk),
      yInitial(Nk),
      x05(Nk),
      y05(Nk),
      u05(Nk),
      v05(Nk),
      e05(Nt),
      Fu(Nk),
      Fv(Nk),
      Mkx(Nk, Nk),
      Mky(Nk, Nk),
      Mt_inv(Nt, Nt),
      Fx(Nk, Nt),
      Fy(Nk, Nt) {
  assert(order > 0);
  initInitializers();
  initBasisValues();
  initKinematicVectors();
  initThermodynamicVector();
  initKinematicMassMatrix();
  initThermodynamicInverseMassMatrix();
  initForceMatrices();
  initSolvers();
}

void FEMALEMethod::calc() {
  RK2step();
  // t += dt;
}
void FEMALEMethod::calcdt() const {}

void FEMALEMethod::dumpData() const {
  std::string filename =
      std::format("{}_{:.3f}.txt", problem.name, t * problem.tMul);
  if (name != problem.name) {
    filename = name + "/" + filename;
  }
  auto outputFunction = [this](std::ofstream ofs) {
    const size_t imax = xSize * order;
    const size_t jmax = ySize * order;
    const double celldx = (problem.xmax - problem.xmin) / xSize;
    const double celldy = (problem.ymax - problem.ymin) / ySize;
    ofs << imax << " " << jmax << std::endl;
    for (size_t i = 0; i < imax; i++) {
      for (size_t j = 0; j < jmax; j++) {
        size_t celli = i / order;
        size_t cellj = j / order;
        double xij = 0.0;
        double yij = 0.0;
        double uij = 0.0;
        double vij = 0.0;

        Eigen::Matrix<double, 2, 2> jacobian =
            Eigen::Matrix<double, 2, 2>::Zero();
        Eigen::Matrix<double, 2, 2> jacobianInitial =
            Eigen::Matrix<double, 2, 2>::Zero();
        for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
          size_t indK = getKinematicIndexFromCell(celli, cellj, k);
          double xk = x[indK];
          double yk = y[indK];
          double uk = u[indK];
          double vk = v[indK];
          double xkInitial = xInitial[indK];
          double ykInitial = yInitial[indK];
          double basis2D =
              output1DKinematicValues[(k % (order + 1)) * order + i % order] *
              output1DKinematicValues[(k / (order + 1)) * order + j % order];
          double basis2Ddx =
              output1DdxKinematicValues[(k % (order + 1)) * order + i % order] *
              output1DKinematicValues[(k / (order + 1)) * order + j % order];
          double basis2Ddy =
              output1DKinematicValues[(k % (order + 1)) * order + i % order] *
              output1DdxKinematicValues[(k / (order + 1)) * order + j % order];
          jacobian(0, 0) += xk * basis2Ddx;
          jacobian(0, 1) += xk * basis2Ddy;
          jacobian(1, 0) += yk * basis2Ddx;
          jacobian(1, 1) += yk * basis2Ddy;
          jacobianInitial(0, 0) += xkInitial * basis2Ddx;
          jacobianInitial(0, 1) += xkInitial * basis2Ddy;
          jacobianInitial(1, 0) += ykInitial * basis2Ddx;
          jacobianInitial(1, 1) += ykInitial * basis2Ddy;
          xij += xk * basis2D;
          yij += yk * basis2D;
          uij += uk * basis2D;
          vij += vk * basis2D;
        }

        double eij = 0.0;
        for (size_t k = 0; k < order * order; k++) {
          size_t indK = getThermodynamicIndexFromCell(celli, cellj, k);
          double ek = e[indK];
          eij += ek * output1DThermoValues[(k % order) * order + i % order] *
                 output1DThermoValues[(k / order) * order + j % order];
        }

        double xlocal = (0.5 + i % order) / order;
        double ylocal = (0.5 + j % order) / order;
        double rhoInitial =
            rhoInitializer(problem.xmin + celldx * (celli + xlocal),
                           problem.ymin + celldy * (cellj + ylocal));
        double rhoij =
            rhoInitial * jacobianInitial.determinant() / jacobian.determinant();
        auto eos = eosInitializer(problem.xmin + celldx * (celli + xlocal),
                                  problem.ymin + celldy * (cellj + ylocal));
        double pij = eos->getp(rhoij, eij);

        std::string line = std::format(
            "{} {} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e} {:.6e}", i, j, xij,
            yij, uij, vij, rhoij, pij, eij);
        ofs << line << std::endl;
      }
    }
  };
  writer.dumpData(outputFunction, filename);
}

void FEMALEMethod::dumpGrid() const {
  std::string filename =
      std::format("{}_grid_{:.3f}.txt", problem.name, t * problem.tMul);
  if (name != problem.name) {
    filename = name + "/" + filename;
  }
  auto outputFunction = [this](std::ofstream ofs) {
    const size_t imax = xSize * order + 1;
    const size_t jmax = ySize * order + 1;
    ofs << imax << " " << jmax << std::endl;
    for (size_t i = 0; i < imax; i++) {
      for (size_t j = 0; j < jmax; j++) {
        double xij = x(i * jmax + j);
        double yij = y(i * jmax + j);

        std::string line = std::format("{} {} {:.6e} {:.6e}", i, j, xij, yij);
        ofs << line << std::endl;
      }
    }
  };
  writer.dumpData(outputFunction, filename);
}

double FEMALEMethod::quadKinematicCellMass(size_t celli, size_t cellj,
                                           size_t basisi, size_t basisj) {
  double output = 0.0;
  const double celldx = (problem.xmax - problem.xmin) / xSize;
  const double celldy = (problem.ymax - problem.ymin) / ySize;
  const size_t indMin = getLegendreStartIndex(kinematicMassQuadOrder);
  for (size_t i = 0; i < kinematicMassQuadOrder; i++) {
    for (size_t j = 0; j < kinematicMassQuadOrder; j++) {
      double xlocal = legendreAbscissas[indMin + i];
      double ylocal = legendreAbscissas[indMin + j];
      double xij = problem.xmin + celldx * (celli + xlocal);
      double yij = problem.ymin + celldy * (cellj + ylocal);

      double rho = rhoInitializer(xij, yij);
      double basisiValue = kinematicMass1DValues[(basisi % (order + 1)) *
                                                     kinematicMassQuadOrder +
                                                 i] *
                           kinematicMass1DValues[(basisi / (order + 1)) *
                                                     kinematicMassQuadOrder +
                                                 j];
      double basisjValue = kinematicMass1DValues[(basisj % (order + 1)) *
                                                     kinematicMassQuadOrder +
                                                 i] *
                           kinematicMass1DValues[(basisj / (order + 1)) *
                                                     kinematicMassQuadOrder +
                                                 j];

      Eigen::Matrix<double, 2, 2> jacobian =
          Eigen::Matrix<double, 2, 2>::Zero();
      for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
        size_t indK = getKinematicIndexFromCell(celli, cellj, k);
        double xk = x(indK);
        double yk = y(indK);
        double basis2Ddx =
            kinematicMass1DdxValues[(k % (order + 1)) * kinematicMassQuadOrder +
                                    i] *
            kinematicMass1DValues[(k / (order + 1)) * kinematicMassQuadOrder +
                                  j];
        double basis2Ddy =
            kinematicMass1DValues[(k % (order + 1)) * kinematicMassQuadOrder +
                                  i] *
            kinematicMass1DdxValues[(k / (order + 1)) * kinematicMassQuadOrder +
                                    j];
        jacobian(0, 0) += xk * basis2Ddx;
        jacobian(0, 1) += xk * basis2Ddy;
        jacobian(1, 0) += yk * basis2Ddx;
        jacobian(1, 1) += yk * basis2Ddy;
      }
      double jacDet = std::abs(jacobian.determinant());

      output += legendreWeights[indMin + i] * legendreWeights[indMin + j] *
                rho * basisiValue * basisjValue * jacDet;
    }
  }
  return output;
}

double FEMALEMethod::quadThermoCellMass(size_t celli, size_t cellj,
                                        size_t basisi, size_t basisj) {
  double output = 0.0;
  const size_t indMin = getLegendreStartIndex(thermoMassQuadOrder);
  const double celldx = (problem.xmax - problem.xmin) / xSize;
  const double celldy = (problem.ymax - problem.ymin) / ySize;
  for (size_t i = 0; i < thermoMassQuadOrder; i++) {
    for (size_t j = 0; j < thermoMassQuadOrder; j++) {
      double xlocal = legendreAbscissas[indMin + i];
      double ylocal = legendreAbscissas[indMin + j];
      double xij = problem.xmin + celldx * (celli + xlocal);
      double yij = problem.ymin + celldy * (cellj + ylocal);

      double rho = rhoInitializer(xij, yij);
      double basisiValue =
          thermoMass1DThermoValues[(basisi % order) * thermoMassQuadOrder + i] *
          thermoMass1DThermoValues[(basisi / order) * thermoMassQuadOrder + j];
      double basisjValue =
          thermoMass1DThermoValues[(basisj % order) * thermoMassQuadOrder + i] *
          thermoMass1DThermoValues[(basisj / order) * thermoMassQuadOrder + j];

      Eigen::Matrix<double, 2, 2> jacobian =
          Eigen::Matrix<double, 2, 2>::Zero();
      for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
        size_t indK = getKinematicIndexFromCell(celli, cellj, k);
        double xk = x(indK);
        double yk = y(indK);
        double basis2Ddx = thermoMass1DdxKinematicValues
                               [(k % (order + 1)) * thermoMassQuadOrder + i] *
                           thermoMass1DKinematicValues[(k / (order + 1)) *
                                                           thermoMassQuadOrder +
                                                       j];
        double basis2Ddy = thermoMass1DKinematicValues[(k % (order + 1)) *
                                                           thermoMassQuadOrder +
                                                       i] *
                           thermoMass1DdxKinematicValues
                               [(k / (order + 1)) * thermoMassQuadOrder + j];
        jacobian(0, 0) += xk * basis2Ddx;
        jacobian(0, 1) += xk * basis2Ddy;
        jacobian(1, 0) += yk * basis2Ddx;
        jacobian(1, 1) += yk * basis2Ddy;
      }
      double jacDet = std::abs(jacobian.determinant());

      output += legendreWeights[indMin + i] * legendreWeights[indMin + j] *
                rho * basisiValue * basisjValue * jacDet;
    }
  }
  return output;
}

std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 2>
FEMALEMethod::quadForceCell(size_t celli, size_t cellj,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
                            const Eigen::Matrix<double, Eigen::Dynamic, 1> &e) {
  std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 2> output{
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(
          (order + 1) * (order + 1), order * order),
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(
          (order + 1) * (order + 1), order * order)};

  const size_t indMin = getLegendreStartIndex(forceQuadOrder);
  for (size_t i = 0; i < forceQuadOrder; i++) {
    for (size_t j = 0; j < forceQuadOrder; j++) {
      Eigen::Matrix<double, 2, 2> jacobian =
          Eigen::Matrix<double, 2, 2>::Zero();
      for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
        size_t indK = getKinematicIndexFromCell(celli, cellj, k);
        double xk = x(indK);
        double yk = y(indK);
        double basis2Ddx =
            force1DdxKinematicValues[(k % (order + 1)) * forceQuadOrder + i] *
            force1DKinematicValues[(k / (order + 1)) * forceQuadOrder + j];
        double basis2Ddy =
            force1DKinematicValues[(k % (order + 1)) * forceQuadOrder + i] *
            force1DdxKinematicValues[(k / (order + 1)) * forceQuadOrder + j];
        jacobian(0, 0) += xk * basis2Ddx;
        jacobian(0, 1) += xk * basis2Ddy;
        jacobian(1, 0) += yk * basis2Ddx;
        jacobian(1, 1) += yk * basis2Ddy;
      }
      Eigen::JacobiSVD<decltype(jacobian)> svd(jacobian);
      double hmin = svd.singularValues().minCoeff() / order;
      double soundSpeed = 0.0;
      double rhoLocal = 0.0;
      double maxViscosityCoeff = 0.0;
      double jacDet = std::abs(jacobian.determinant());

      Eigen::Matrix<double, 2, 2> stressTensor =
          calcStressTensor(u, v, e, jacobian, soundSpeed, rhoLocal,
                           maxViscosityCoeff, celli, cellj, i, j);

      for (size_t basisKinematic = 0;
           basisKinematic < (order + 1) * (order + 1); basisKinematic++) {
        for (size_t basisThermo = 0; basisThermo < order * order;
             basisThermo++) {
          double basisKinematic2Ddx =
              force1DdxKinematicValues[(basisKinematic % (order + 1)) *
                                           forceQuadOrder +
                                       i] *
              force1DKinematicValues[(basisKinematic / (order + 1)) *
                                         forceQuadOrder +
                                     j];
          double basisKinematic2Ddy =
              force1DKinematicValues[(basisKinematic % (order + 1)) *
                                         forceQuadOrder +
                                     i] *
              force1DdxKinematicValues[(basisKinematic / (order + 1)) *
                                           forceQuadOrder +
                                       j];
          Eigen::Matrix<double, 2, 2> gradBasisx{
              {basisKinematic2Ddx, basisKinematic2Ddy}, {0.0, 0.0}};
          Eigen::Matrix<double, 2, 2> gradBasisy{
              {0.0, 0.0}, {basisKinematic2Ddx, basisKinematic2Ddy}};
          double thermobasis =
              force1DThermoValues[(basisThermo % order) * forceQuadOrder + i] *
              force1DThermoValues[(basisThermo / order) * forceQuadOrder + j];

          Eigen::Matrix<double, 2, 2> rhsx = jacobian.inverse();
          Eigen::Matrix<double, 2, 2> rhsy = rhsx;
          rhsx = rhsx * gradBasisx;
          rhsx *= thermobasis * jacDet;
          rhsy = rhsy * gradBasisy;
          rhsy *= thermobasis * jacDet;

          // output[0] += legendreWeights[indMin + i] *
          //              legendreWeights[indMin + j] *
          //              (stressTensor(0, 0) * rhsx(0, 0) +
          //               stressTensor(0, 1) * rhsx(0, 1) +
          //               stressTensor(1, 0) * rhsx(1, 0) +
          //               stressTensor(1, 1) * rhsx(1, 1));
          // output[1] += legendreWeights[indMin + i] *
          //              legendreWeights[indMin + j] *
          //              (stressTensor(0, 0) * rhsy(0, 0) +
          //               stressTensor(0, 1) * rhsy(0, 1) +
          //               stressTensor(1, 0) * rhsy(1, 0) +
          //               stressTensor(1, 1) * rhsy(1, 1));
          output[0](basisKinematic, basisThermo) +=
              legendreWeights[indMin + i] * legendreWeights[indMin + j] *
              (stressTensor.transpose() * rhsx).trace();
          output[1](basisKinematic, basisThermo) +=
              legendreWeights[indMin + i] * legendreWeights[indMin + j] *
              (stressTensor.transpose() * rhsy).trace();
        }
      }
      calcTau(hmin, soundSpeed, rhoLocal, maxViscosityCoeff);
    }
  }
  return output;
}

size_t FEMALEMethod::getKinematicIndexFromCell(const size_t celli,
                                               const size_t cellj,
                                               const size_t k) const {
  assert(k < (order + 1) * (order + 1));
  return (celli * order + k % (order + 1)) * (ySize * order + 1) +
         cellj * order + k / (order + 1);
}
size_t FEMALEMethod::getThermodynamicIndexFromCell(const size_t celli,
                                                   const size_t cellj,
                                                   const size_t k) const {
  assert(k < order * order);
  return celli * order * order * ySize + cellj * order * order + k;
}

void FEMALEMethod::initInitializers() {
  rhoInitializer =
      [xmin = problem.xmin, xmax = problem.xmax, ymin = problem.ymin,
       ymax = problem.ymax, xSize = this->xSize, ySize = this->ySize,
       rhoInitializer = problem.rhoInitializer](double x, double y) {
        double dx = (xmax - xmin) / xSize;
        double dy = (ymax - ymin) / ySize;
        double xij = xmin + (std::floor((x - xmin) / dx) + 0.5) * dx;
        double yij = ymin + (std::floor((y - ymin) / dy) + 0.5) * dy;
        return rhoInitializer(xij, yij);
      };
  eosInitializer =
      [xmin = problem.xmin, xmax = problem.xmax, ymin = problem.ymin,
       ymax = problem.ymax, xSize = this->xSize, ySize = this->ySize,
       eosInitializer = problem.eosInitializer](double x, double y) {
        double dx = (xmax - xmin) / xSize;
        double dy = (ymax - ymin) / ySize;
        double xij = xmin + (std::floor((x - xmin) / dx) + 0.5) * dx;
        double yij = ymin + (std::floor((y - ymin) / dy) + 0.5) * dy;
        return eosInitializer(xij, yij);
      };
}

void FEMALEMethod::initBasisValues() {
  initKinematicMassBasisValues();
  initThermodynamicMassBasisValues();
  initForceBasisValues();
  initOutputBasisValues();
}

void FEMALEMethod::initKinematicMassBasisValues() {
  kinematicMass1DValues.reserve((order + 1) * kinematicMassQuadOrder);
  kinematicMass1DdxValues.reserve((order + 1) * kinematicMassQuadOrder);

  const size_t indMin = getLegendreStartIndex(kinematicMassQuadOrder);
  for (size_t k = 0; k < order + 1; k++) {
    for (size_t i = indMin; i < indMin + kinematicMassQuadOrder; i++) {
      double xi = legendreAbscissas[i];
      kinematicMass1DValues.push_back(lobattoBasis1D(xi, order, k));
      kinematicMass1DdxValues.push_back(lobattoBasis1Ddx(xi, order, k));
    }
  }
}
void FEMALEMethod::initThermodynamicMassBasisValues() {
  thermoMass1DKinematicValues.reserve((order + 1) * thermoMassQuadOrder);
  thermoMass1DdxKinematicValues.reserve((order + 1) * thermoMassQuadOrder);
  thermoMass1DThermoValues.reserve(order * thermoMassQuadOrder);

  const size_t indMin = getLegendreStartIndex(thermoMassQuadOrder);
  for (size_t k = 0; k < order + 1; k++) {
    for (size_t i = indMin; i < indMin + thermoMassQuadOrder; i++) {
      double xi = legendreAbscissas[i];
      thermoMass1DKinematicValues.push_back(lobattoBasis1D(xi, order, k));
      thermoMass1DdxKinematicValues.push_back(lobattoBasis1Ddx(xi, order, k));
    }
  }

  for (size_t k = 0; k < order; k++) {
    for (size_t i = indMin; i < indMin + thermoMassQuadOrder; i++) {
      double xi = legendreAbscissas[i];
      thermoMass1DThermoValues.push_back(lobattoBasis1D(xi, order - 1, k));
    }
  }
}
void FEMALEMethod::initForceBasisValues() {
  force1DKinematicValues.reserve((order + 1) * forceQuadOrder);
  force1DdxKinematicValues.reserve((order + 1) * forceQuadOrder);

  force1DThermoValues.reserve(order * forceQuadOrder);

  const size_t indMin = getLegendreStartIndex(forceQuadOrder);
  for (size_t k = 0; k < order + 1; k++) {
    for (size_t i = indMin; i < indMin + forceQuadOrder; i++) {
      double xi = legendreAbscissas[i];
      force1DKinematicValues.push_back(lobattoBasis1D(xi, order, k));
      force1DdxKinematicValues.push_back(lobattoBasis1Ddx(xi, order, k));
    }
  }

  for (size_t k = 0; k < order; k++) {
    for (size_t i = indMin; i < indMin + forceQuadOrder; i++) {
      double xi = legendreAbscissas[i];
      force1DThermoValues.push_back(lobattoBasis1D(xi, order - 1, k));
    }
  }
}

void FEMALEMethod::initOutputBasisValues() {
  output1DKinematicValues.reserve((order + 1) * order);
  output1DdxKinematicValues.reserve((order + 1) * order);
  output1DThermoValues.reserve(order * order);

  double dx = 1.0 / order;
  for (size_t k = 0; k < order + 1; k++) {
    for (size_t i = 0; i < order; i++) {
      double xi = 0.5 * dx + dx * i;
      output1DKinematicValues.push_back(lobattoBasis1D(xi, order, k));
      output1DdxKinematicValues.push_back(lobattoBasis1Ddx(xi, order, k));
    }
  }

  for (size_t k = 0; k < order; k++) {
    for (size_t i = 0; i < order; i++) {
      double xi = 0.5 * dx + dx * i;
      output1DThermoValues.push_back(lobattoBasis1D(xi, order - 1, k));
    }
  }
}

void FEMALEMethod::initKinematicVectors() {
  const size_t imax = xSize * order + 1;  // number of nodes in x direction
  const size_t jmax = ySize * order + 1;  // number of nodes in y direction
  const double celldx = (problem.xmax - problem.xmin) / xSize;
  const double celldy = (problem.ymax - problem.ymin) / ySize;
  switch (problem.dimension) {
    case ProblemDimension::e1D:
      l0 = celldx / order;
      break;
    case ProblemDimension::e2D:
      l0 = std::sqrt(celldx * celldy) / order;
      break;
  }
  for (size_t i = 0; i < imax; i++) {
    for (size_t j = 0; j < jmax; j++) {
      size_t indMin = getLobattoStartIndex(order);
      size_t lobattoxIndex = indMin + i % order;
      size_t lobattoyIndex = indMin + j % order;
      double xij = problem.xmin + celldx * (static_cast<double>(i / order) +
                                        lobattoAbscissas[lobattoxIndex]);
      double yij = problem.ymin + celldy * (static_cast<double>(j / order) +
                                        lobattoAbscissas[lobattoyIndex]);
      x(i * jmax + j) = xij;
      xInitial(i * jmax + j) = xij;
      y(i * jmax + j) = yij;
      yInitial(i * jmax + j) = yij;
      u(i * jmax + j) = problem.uInitializer(xij, yij);
      v(i * jmax + j) = problem.vInitializer(xij, yij);
    }
  }
}

void FEMALEMethod::initThermodynamicVector() {
  double celldx = (problem.xmax - problem.xmin) / xSize;
  double celldy = (problem.ymax - problem.ymin) / ySize;
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      double localdx = 0.5;
      double localdy = 0.5;
      for (size_t k = 0; k < order * order; k++) {
        // if (order != 1) {
        //   size_t indMin = getLobattoStartIndex(order - 1);
        //   size_t lobattoxIndex = indMin + k % order;
        //   size_t lobattoyIndex = indMin + k / order;
        //   localdx = lobattoAbscissas[lobattoxIndex];
        //   localdy = lobattoAbscissas[lobattoyIndex];
        // }
        double xij = problem.xmin + celldx * (celli + localdx);
        double yij = problem.ymin + celldy * (cellj + localdy);
        auto eos = eosInitializer(xij, yij);
        auto p = problem.pInitializer(xij, yij);
        auto rho = rhoInitializer(xij, yij);
        e(getThermodynamicIndexFromCell(celli, cellj, k)) = eos->gete(rho, p);
      }
    }
  }
}

void FEMALEMethod::initKinematicMassMatrix() {
  Mkx.reserve(Eigen::VectorXi::Constant(
      Nk, (2 * order + 1) *
              (2 * order + 1)));  // reserve maximum possible space for each row
  Mky.reserve(Eigen::VectorXi::Constant(
      Nk, (2 * order + 1) *
              (2 * order + 1)));  // reserve maximum possible space for each row
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        for (size_t j = i; j < (order + 1) * (order + 1); j++) {
          double quadResult = quadKinematicCellMass(celli, cellj, i, j);

          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
          Mkx.coeffRef(indI, indJ) += quadResult;
          if (indI != indJ) {
            Mkx.coeffRef(indJ, indI) += quadResult;
          }
          Mky.coeffRef(indI, indJ) += quadResult;
          if (indI != indJ) {
            Mky.coeffRef(indJ, indI) += quadResult;
          }
        }
      }
    }
  }
  resolveBoundaryMass();
  Mkx.makeCompressed();
  Mky.makeCompressed();
}

void FEMALEMethod::initThermodynamicInverseMassMatrix() {
  Mt_inv.reserve(Eigen::VectorXi::Constant(Nt, order * order));
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellMt(
          order * order, order * order);
      for (size_t i = 0; i < order * order; i++) {
        for (size_t j = i; j < order * order; j++) {
          double quadResult = quadThermoCellMass(celli, cellj, i, j);
          cellMt(i, j) = quadResult;
          if (i != j) {
            cellMt(j, i) = quadResult;
          }
        }
      }
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellMt_inv(
          order * order, order * order);
      cellMt_inv = cellMt.inverse();
      for (size_t i = 0; i < order * order; i++) {
        for (size_t j = 0; j < order * order; j++) {
          size_t indI = getThermodynamicIndexFromCell(celli, cellj, i);
          size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
          Mt_inv.insert(indI, indJ) = cellMt_inv(i, j);
        }
      }
    }
  }
  Mt_inv.makeCompressed();
}

void FEMALEMethod::initForceMatrices() {
  Fx.reserve(Eigen::VectorXi::Constant(Nt, (order + 1) * (order + 1)));
  Fy.reserve(Eigen::VectorXi::Constant(Nt, (order + 1) * (order + 1)));
}

void FEMALEMethod::initSolvers() {
  kinematicSolverx.compute(Mkx);
  assert(kinematicSolverx.info() == Eigen::Success);
  kinematicSolvery.compute(Mky);
  assert(kinematicSolvery.info() == Eigen::Success);
}

void FEMALEMethod::calcForceMatrices(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &e) {
#pragma omp parallel for
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      std::array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, 2>
          quadResult = quadForceCell(celli, cellj, x, y, u, v, e);
      for (size_t kinematikk = 0; kinematikk < (order + 1) * (order + 1);
           kinematikk++) {
        for (size_t thermok = 0; thermok < order * order; thermok++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, kinematikk);
          size_t indJ = getThermodynamicIndexFromCell(celli, cellj, thermok);

          Fx.coeffRef(indI, indJ) = quadResult[0](kinematikk, thermok);
          Fy.coeffRef(indI, indJ) = quadResult[1](kinematikk, thermok);
        }
      }
    }
  }
  resolveBoundaryForce();
}

Eigen::Matrix<double, 2, 2> FEMALEMethod::calcStressTensor(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &e,
    const Eigen::Matrix<double, 2, 2> &jacobian, double &soundSpeed,
    double &rhoLocal, double &maxViscosityCoeff, const size_t celli,
    const size_t cellj, const size_t i, const size_t j) const {
  Eigen::Matrix<double, 2, 2> stressTensor =
      Eigen::Matrix<double, 2, 2>::Zero();

  const size_t indMin = getLegendreStartIndex(forceQuadOrder);
  const double celldx = (problem.xmax - problem.xmin) / xSize;
  const double celldy = (problem.ymax - problem.ymin) / ySize;
  const double xlocal = legendreAbscissas[indMin + i];
  const double ylocal = legendreAbscissas[indMin + j];
  const double xij = problem.xmin + celldx * (celli + xlocal);
  const double yij = problem.ymin + celldy * (cellj + ylocal);

  double eLocal = 0.0;
  for (size_t k = 0; k < order * order; k++) {
    size_t indK = getThermodynamicIndexFromCell(celli, cellj, k);
    eLocal += e[indK] * force1DThermoValues[(k % order) * forceQuadOrder + i] *
              force1DThermoValues[(k / order) * forceQuadOrder + j];
  }

  Eigen::Matrix<double, 2, 2> jacobianInitial =
      Eigen::Matrix<double, 2, 2>::Zero();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    size_t indK = getKinematicIndexFromCell(celli, cellj, k);
    double xkInitial = xInitial(indK);
    double ykInitial = yInitial(indK);
    double basis2Ddx =
        force1DdxKinematicValues[(k % (order + 1)) * forceQuadOrder + i] *
        force1DKinematicValues[(k / (order + 1)) * forceQuadOrder + j];
    double basis2Ddy =
        force1DKinematicValues[(k % (order + 1)) * forceQuadOrder + i] *
        force1DdxKinematicValues[(k / (order + 1)) * forceQuadOrder + j];
    jacobianInitial(0, 0) += xkInitial * basis2Ddx;
    jacobianInitial(0, 1) += xkInitial * basis2Ddy;
    jacobianInitial(1, 0) += ykInitial * basis2Ddx;
    jacobianInitial(1, 1) += ykInitial * basis2Ddy;
  }

  double rhoInitial = rhoInitializer(xij, yij);
  rhoLocal = rhoInitial *
             std::abs(jacobianInitial.determinant() / jacobian.determinant());

  auto eos = eosInitializer(xij, yij);
  double pLocal = eos->getp(rhoLocal, eLocal);
  soundSpeed = eos->getc(rhoLocal, pLocal);

  stressTensor = -pLocal * Eigen::Matrix<double, 2, 2>::Identity();
  stressTensor +=
      calcArtificialViscosity(u, v, jacobian, jacobianInitial, soundSpeed,
                              rhoLocal, maxViscosityCoeff, celli, cellj, i, j);

  return stressTensor;
}

Eigen::Matrix<double, 2, 2> FEMALEMethod::calcArtificialViscosity(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &u,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &v,
    const Eigen::Matrix<double, 2, 2> &jacobian,
    const Eigen::Matrix<double, 2, 2> &jacobianInitial, double soundSpeed,
    double rhoLocal, double &maxViscosityCoeff, const size_t celli,
    const size_t cellj, const size_t i, const size_t j) const {
  Eigen::Matrix<double, 2, 2> output = Eigen::Matrix<double, 2, 2>::Zero();

  Eigen::Matrix<double, 2, 2> velocityGrad =
      Eigen::Matrix<double, 2, 2>::Zero();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    size_t indK = getKinematicIndexFromCell(celli, cellj, k);
    double uk = u[indK];
    double vk = v[indK];
    double basis2Ddx =
        force1DdxKinematicValues[(k % (order + 1)) * forceQuadOrder + i] *
        force1DKinematicValues[(k / (order + 1)) * forceQuadOrder + j];
    double basis2Ddy =
        force1DKinematicValues[(k % (order + 1)) * forceQuadOrder + i] *
        force1DdxKinematicValues[(k / (order + 1)) * forceQuadOrder + j];
    velocityGrad(0, 0) += uk * basis2Ddx;
    velocityGrad(0, 1) += uk * basis2Ddy;
    velocityGrad(1, 0) += vk * basis2Ddx;
    velocityGrad(1, 1) += vk * basis2Ddy;
  }
  velocityGrad = jacobian.inverse() * velocityGrad;
  Eigen::Matrix<double, 2, 2> symVelocityGrad =
      0.5 * (velocityGrad + velocityGrad.transpose());

  Eigen::SelfAdjointEigenSolver<decltype(symVelocityGrad)> eigensolver(
      symVelocityGrad);

  double e1 = eigensolver.eigenvalues()(0);
  double e2 = eigensolver.eigenvalues()(1);
  Eigen::Matrix<double, 2, 1> v1 = eigensolver.eigenvectors().col(0);
  Eigen::Matrix<double, 2, 1> v2 = eigensolver.eigenvectors().col(1);
  double viscosityCoeff1 = calcViscosityCoeff(
      jacobian, jacobianInitial, velocityGrad, e1, v1, soundSpeed, rhoLocal);
  double viscosityCoeff2 = calcViscosityCoeff(
      jacobian, jacobianInitial, velocityGrad, e2, v2, soundSpeed, rhoLocal);

  Eigen::Matrix<double, 2, 2> v1Tensor{{v1(0) * v1(0), v1(0) * v1(1)},
                                       {v1(0) * v1(1), v1(1) * v1(1)}};
  Eigen::Matrix<double, 2, 2> v2Tensor{{v2(0) * v2(0), v2(0) * v2(1)},
                                       {v2(0) * v2(1), v2(1) * v2(1)}};
  output += viscosityCoeff1 * e1 * v1Tensor;
  output += viscosityCoeff2 * e2 * v2Tensor;

  maxViscosityCoeff = std::max(viscosityCoeff1, viscosityCoeff2);

  return output;
}

double FEMALEMethod::calcViscosityCoeff(
    const Eigen::Matrix<double, 2, 2> &jacobian,
    const Eigen::Matrix<double, 2, 2> &jacobianInitial,
    const Eigen::Matrix<double, 2, 2> &velocityGrad, double eigenvalue,
    const Eigen::Matrix<double, 2, 1> &eigenvector, double soundSpeed,
    double rhoLocal) const {
  double psi0 = 0.0;
  double velgradNorm = velocityGrad.norm();
  if (velgradNorm != 0.0) {
    double velocityScalarGrad = velocityGrad(0, 0) + velocityGrad(1, 1);
    velocityScalarGrad = std::abs(velocityScalarGrad);

    psi0 = velocityScalarGrad / velgradNorm;
  }
  double psi1 = eigenvalue < 0 ? 1 : 0;

  Eigen::Matrix<double, 2, 2> jacobianInitial_inv = jacobianInitial.inverse();
  Eigen::Matrix<double, 2, 2> jacobianMapping = jacobianInitial_inv * jacobian;
  double localLengthScale =
      l0 * (jacobianMapping * eigenvector).norm() / eigenvector.norm();

  double lhs = q2 * localLengthScale * localLengthScale * std::abs(eigenvalue);
  double rhs = q1 * psi0 * psi1 * localLengthScale * soundSpeed;
  double output = rhoLocal * (lhs + rhs);
  return output;
}

void FEMALEMethod::calcTau(double hmin, double soundSpeed, double rhoLocal,
                           double maxViscosityCoeff) {
  double denominator = soundSpeed / hmin +
                       alphamu * maxViscosityCoeff / (rhoLocal * hmin * hmin);
  double tauLocal = alpha / denominator;
#pragma omp critical(tauUpdate)
  {
    if (tauLocal < tau) {
      tau = tauLocal;
    }
  }
}

void FEMALEMethod::resolveBoundaryMass() {
  resolveLeftBoundaryMass();
  resolveTopBoundaryMass();
  resolveRightBoundaryMass();
  resolveBottomBoundaryMass();
}
void FEMALEMethod::resolveLeftBoundaryMass() {
  const size_t celli = 0;
  switch (problem.leftBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = 0; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mkx.coeffRef(indI, indJ) = 0.0;
            Mkx.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = 0; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mkx.coeffRef(indI, indJ) = 0.0;
            Mkx.coeffRef(indJ, indI) = 0.0;
            Mky.coeffRef(indI, indJ) = 0.0;
            Mky.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveTopBoundaryMass() {
  const size_t cellj = ySize - 1;
  switch (problem.topBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = order * (order + 1); i < (order + 1) * (order + 1);
             i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mky.coeffRef(indI, indJ) = 0.0;
            Mky.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = order * (order + 1); i < (order + 1) * (order + 1);
             i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mkx.coeffRef(indI, indJ) = 0.0;
            Mkx.coeffRef(indJ, indI) = 0.0;
            Mky.coeffRef(indI, indJ) = 0.0;
            Mky.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveRightBoundaryMass() {
  const size_t celli = xSize - 1;
  switch (problem.rightBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = order; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mkx.coeffRef(indI, indJ) = 0.0;
            Mkx.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = order; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mkx.coeffRef(indI, indJ) = 0.0;
            Mkx.coeffRef(indJ, indI) = 0.0;
            Mky.coeffRef(indI, indJ) = 0.0;
            Mky.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveBottomBoundaryMass() {
  const size_t cellj = 0;
  switch (problem.bottomBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = 0; i < order + 1; i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mky.coeffRef(indI, indJ) = 0.0;
            Mky.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = 0; i < order + 1; i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            Mkx.coeffRef(indI, indJ) = 0.0;
            Mkx.coeffRef(indJ, indI) = 0.0;
            Mky.coeffRef(indI, indJ) = 0.0;
            Mky.coeffRef(indJ, indI) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveBoundaryForce() {
  resolveLeftBoundaryForce();
  resolveTopBoundaryForce();
  resolveRightBoundaryForce();
  resolveBottomBoundaryForce();
}
void FEMALEMethod::resolveLeftBoundaryForce() {
  const size_t celli = 0;
  switch (problem.leftBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = 0; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fx.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = 0; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fx.coeffRef(indI, indJ) = 0.0;
            Fy.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveTopBoundaryForce() {
  const size_t cellj = ySize - 1;
  switch (problem.topBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = order * (order + 1); i < (order + 1) * (order + 1);
             i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fy.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = order * (order + 1); i < (order + 1) * (order + 1);
             i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fx.coeffRef(indI, indJ) = 0.0;
            Fy.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveRightBoundaryForce() {
  const size_t celli = xSize - 1;
  switch (problem.rightBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = order; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fx.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = order; i < (order + 1) * (order + 1); i += order + 1) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fx.coeffRef(indI, indJ) = 0.0;
            Fy.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}
void FEMALEMethod::resolveBottomBoundaryForce() {
  const size_t cellj = 0;
  switch (problem.bottomBoundaryType) {
    case BoundaryType::eFree:
      break;
    case BoundaryType::eWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = 0; i < order + 1; i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fy.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eNoSlipWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = 0; i < order + 1; i++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < order * order; j++) {
            size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
            Fx.coeffRef(indI, indJ) = 0.0;
            Fy.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}

void FEMALEMethod::RK2step() {
  // dt / 2
  while (true) {
    tau = std::numeric_limits<double>::max();
    calcForceMatrices(x, y, u, v, e);
    if (dt >= tau) {
      dt = beta1 * tau;
    }

    Fu = Fx * Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(Nt);
    Fv = Fy * Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(Nt);

    u05 = kinematicSolverx.solve(Fu);
    assert(kinematicSolverx.info() == Eigen::Success);
    v05 = kinematicSolvery.solve(Fv);
    assert(kinematicSolvery.info() == Eigen::Success);
    u05 *= -0.5 * dt;
    u05 += u;
    v05 *= -0.5 * dt;
    v05 += v;

    e05 = Mt_inv * (Fx.transpose() * u05 + Fy.transpose() * v05);
    e05 *= 0.5 * dt;
    e05 += e;

    x05 = x + 0.5 * dt * u05;
    y05 = y + 0.5 * dt * v05;

    // dt
    calcForceMatrices(x05, y05, u05, v05, e05);
    if (dt >= tau) {
      dt = beta1 * tau;
      continue;
    } else {
      break;
    }
  }

  Fu = Fx * Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(Nt);
  Fv = Fy * Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(Nt);

  u05 = kinematicSolverx.solve(Fu);
  assert(kinematicSolverx.info() == Eigen::Success);
  v05 = kinematicSolvery.solve(Fv);
  assert(kinematicSolvery.info() == Eigen::Success);
  u05 *= -dt;
  u05 += u;
  u.swap(u05);
  u05 += u;
  u05 *= 0.5;

  v05 *= -dt;
  v05 += v;
  v.swap(v05);
  v05 += v;
  v05 *= 0.5;

  e05 = Mt_inv * (Fx.transpose() * u05 + Fy.transpose() * v05);
  e05 *= dt;
  e05 += e;

  x05 = x + dt * u05;
  y05 = y + dt * v05;

  t += dt;
  if (dt <= gamma * tau) {
    dt = beta2 * dt;
  }

  x = x05;
  y = y05;
  e = e05;
}
