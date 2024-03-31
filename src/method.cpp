#include "method.hpp"

#include <algorithm>
#include <cassert>
#include <format>

#include "nodes.hpp"

FEMALEMethod::FEMALEMethod(const std::string &name, const Problem &problem_,
                           size_t xSize_, size_t ySize_, size_t order_)
    : Method(name, problem_),
      xSize(xSize_),
      ySize(ySize_),
      order(order_),
      quadOrder(order * 2),
      Nk((order * xSize + 1) * (order * ySize + 1)),
      Nt((order * xSize) * (order * ySize)),
      Na((order + 1) * xSize * (order + 1) * ySize),
      x(Nk),
      y(Nk),
      u(Nk),
      v(Nk),
      e(Nt),
      rhoQuad((quadOrder * xSize) * (quadOrder * ySize)),
      xInitial(Nk),
      yInitial(Nk),
      x05(Nk),
      y05(Nk),
      u05(Nk),
      v05(Nk),
      e05(Nt),
      Fu(Nk),
      Fv(Nk),
      Fe(Nt),
      xOptimal(Nk),
      yOptimal(Nk),
      uMesh(Nk),
      vMesh(Nk),
      rhoRemap(Na),
      rhoRemapHigh(Na),
      rhoRemapLow(Na),
      rhoRemapAvg(Na),
      rhoRemapMin(Na),
      rhoRemapMax(Na),
      rhoERemap(Na),
      rhoERemapHigh(Na),
      rhoERemapLow(Na),
      rhoERemapAvg(Na),
      eRemapMin(Na),
      eRemapMax(Na),
      Mkx(Nk, Nk),
      Mky(Nk, Nk),
      Mt_inv(Nt, Nt),
      Fx(Nk, Nt),
      Fy(Nk, Nt),
      vectorLaplacianX(Nk, Nk),
      vectorLaplacianY(Nk, Nk),
      vectorMass(Nk, Nk),
      Mv(Nk, Nk),
      M(Na, Na),
      Kv(Nk, Nk),
      K(Na, Na),
      Mlumped(Na, Na),
      MlumpedInv(Na, Na),
      MInv(Na, Na),
      L(Na, Na),
      Kupwinded(Na, Na),
      D(Na, Na),
      antidiffusiveFlux(Na, Na),
      fluxLimitingFactors(Na, Na)

{
  assert(order > 0);
  initInitializers();
  initBasisValues();
  initRhoValues();
  initKinematicVectors();
  initThermodynamicVector();
  initKinematicMassMatrix();
  initThermodynamicInverseMassMatrix();
  initForceMatrices();
  initVectorMatrices();
  initRemapMatrices();
  initSolvers();
}

void FEMALEMethod::calc() {
  RK2step();
  if (++iteration % remapFrequency == 0) {
    remap();
  }
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
        // double rhoInitial =
        //     rhoInitializer(problem.xmin + celldx * (celli + xlocal),
        //                    problem.ymin + celldy * (cellj + ylocal));
        double rho = 0.0;
        for (size_t k = 0; k < quadOrder * quadOrder; k++) {
          const size_t indRho = getQuadIndexFromCell(celli, cellj, k);
          const double rhok = rhoQuad(indRho);
          const double basis2D =
              output1DRhoValues[(k % quadOrder) * order + i % order] *
              output1DRhoValues[(k / quadOrder) * order + j % order];

          rho += rhok * basis2D;
        }
        double rhoij =
            rho * jacobianInitial.determinant() / jacobian.determinant();
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

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
FEMALEMethod::quadKinematicCellMass(size_t celli, size_t cellj) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output(
      (order + 1) * (order + 1), (order + 1) * (order + 1));

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const size_t quadi = quad % quadOrder;
    const size_t quadj = quad / quadOrder;

    const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
    const double rhoLocal = rhoQuad(indRho);

    Eigen::Matrix<double, 2, 2> jacobian =
        getCellJacobian(celli, cellj, quadi, quadj, x, y);
    const double jacDet = std::abs(jacobian.determinant());

    for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
      const double basisi2D =
          kinematicBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                     quadi] *
          kinematicBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                     quadj];
      for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
        const double basisj2D =
            kinematicBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                       quadi] *
            kinematicBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                       quadj];

        output(basisi, basisj) = legendreWeights[indMin + quadi] *
                                 legendreWeights[indMin + quadj] * rhoLocal *
                                 basisi2D * basisj2D * jacDet;
      }
    }
  }
  return output;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
FEMALEMethod::quadThermoCellMass(size_t celli, size_t cellj) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output(order * order,
                                                               order * order);

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const size_t quadi = quad % quadOrder;
    const size_t quadj = quad / quadOrder;

    Eigen::Matrix<double, 2, 2> jacobian =
        getCellJacobian(celli, cellj, quadi, quadj, x, y);
    const double jacDet = std::abs(jacobian.determinant());

    const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
    const double rhoLocal = rhoQuad(indRho);

    for (size_t basisi = 0; basisi < order * order; basisi++) {
      const double basisi2D =
          thermoBasis1DQuadValues[(basisi % order) * quadOrder + quadi] *
          thermoBasis1DQuadValues[(basisi / order) * quadOrder + quadj];
      for (size_t basisj = 0; basisj < order * order; basisj++) {
        const double basisj2D =
            thermoBasis1DQuadValues[(basisj % order) * quadOrder + quadi] *
            thermoBasis1DQuadValues[(basisj / order) * quadOrder + quadj];

        output(basisi, basisj) += legendreWeights[indMin + quadi] *
                                  legendreWeights[indMin + quadj] * rhoLocal *
                                  basisi2D * basisj2D * jacDet;
      }
    }
  }
  return output;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
FEMALEMethod::quadCellVectorMass() {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output{
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(
          (order + 1) * (order + 1), (order + 1) * (order + 1))};
  const size_t indMin = getLegendreStartIndex(quadOrder);

  for (size_t i = 0; i < quadOrder; i++) {
    for (size_t j = 0; j < quadOrder; j++) {
      for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
        for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
          double basisiValue =
              kinematicBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                         i] *
              kinematicBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                         j];
          double basisjValue =
              kinematicBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                         i] *
              kinematicBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                         j];

          output(basisi, basisj) += legendreWeights[indMin + i] *
                                    legendreWeights[indMin + j] * basisiValue *
                                    basisjValue;
        }
      }
    }
  }

  return output;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
FEMALEMethod::quadCellVectorLaplacian() {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output{
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(
          (order + 1) * (order + 1), (order + 1) * (order + 1))};
  const size_t indMin = getLegendreStartIndex(quadOrder);

  for (size_t i = 0; i < quadOrder; i++) {
    for (size_t j = 0; j < quadOrder; j++) {
      for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
        double basisi2Ddx =
            kinematicBasis1DdxQuadValues[(basisi % (order + 1)) * quadOrder +
                                         i] *
            kinematicBasis1DQuadValues[(basisi / (order + 1)) * quadOrder + j];
        double basisi2Ddy =
            kinematicBasis1DQuadValues[(basisi % (order + 1)) * quadOrder + i] *
            kinematicBasis1DdxQuadValues[(basisi / (order + 1)) * quadOrder +
                                         j];
        Eigen::Matrix<double, 2, 2> gradBasisi{{basisi2Ddx, basisi2Ddy},
                                               {basisi2Ddx, basisi2Ddy}};

        for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
          double basisj2Ddx =
              kinematicBasis1DdxQuadValues[(basisj % (order + 1)) * quadOrder +
                                           i] *
              kinematicBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                         j];
          double basisj2Ddy =
              kinematicBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                         i] *
              kinematicBasis1DdxQuadValues[(basisj / (order + 1)) * quadOrder +
                                           j];
          Eigen::Matrix<double, 2, 2> gradBasisj{{basisj2Ddx, basisj2Ddy},
                                                 {basisj2Ddx, basisj2Ddy}};

          double value = (gradBasisi.transpose() * gradBasisj).trace();

          output(basisi, basisj) +=
              legendreWeights[indMin + i] * legendreWeights[indMin + j] * value;
        }
      }
    }
  }

  return output;
}
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> FEMALEMethod::quadCellMv(
    size_t celli, size_t cellj,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output(
      (order + 1) * (order + 1), (order + 1) * (order + 1));

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const size_t quadi = quad % quadOrder;
    const size_t quadj = quad / quadOrder;
    Eigen::Matrix<double, 2, 2> jacobian =
        getCellJacobian(celli, cellj, quadi, quadj, x, y);
    double jacDet = std::abs(jacobian.determinant());

    double rhoLocal = 0.0;
    for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
      const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
      const size_t rhoK = rhoRemap(indK);
      const double basisk2D =
          advectionBasis1DQuadValues[(k % (order + 1)) * quadOrder + quadi] *
          advectionBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];

      rhoLocal += rhoK * basisk2D;
    }

    for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
      double basisiValue =
          kinematicBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                     quadi] *
          kinematicBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                     quadj];
      for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
        double basisjValue =
            kinematicBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                       quadi] *
            kinematicBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                       quadj];

        output(basisi, basisj) += legendreWeights[indMin + quadi] *
                                  legendreWeights[indMin + quadj] * rhoLocal *
                                  basisiValue * basisjValue * jacDet;
      }
    }
  }

  return output;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> FEMALEMethod::quadCellM(
    size_t celli, size_t cellj,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output(
      (order + 1) * (order + 1), (order + 1) * (order + 1));

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quadi = 0; quadi < quadOrder; quadi++) {
    for (size_t quadj = 0; quadj < quadOrder; quadj++) {
      Eigen::Matrix<double, 2, 2> jacobian =
          getCellJacobian(celli, cellj, quadi, quadj, x, y);
      double jacDet = std::abs(jacobian.determinant());

      for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
        double basisi2D =
            advectionBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                       quadi] *
            advectionBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                       quadj];
        for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
          double basisj2D =
              advectionBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                         quadi] *
              advectionBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                         quadj];

          output(basisi, basisj) += legendreWeights[indMin + quadi] *
                                    legendreWeights[indMin + quadj] * basisi2D *
                                    basisj2D * jacDet;
        }
      }
    }
  }
  return output;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> FEMALEMethod::quadCellKv(
    size_t celli, size_t cellj,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output(
      (order + 1) * (order + 1), (order + 1) * (order + 1));

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quadi = 0; quadi < quadOrder; quadi++) {
    for (size_t quadj = 0; quadj < quadOrder; quadj++) {
      Eigen::Matrix<double, 2, 2> jacobian =
          getCellJacobian(celli, cellj, quadi, quadj, x, y);
      const double jacDet = std::abs(jacobian.determinant());

      double rhoLocal = 0.0;
      for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
        const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
        const size_t rhoK = rhoRemap(indK);
        const size_t basisk2D =
            advectionBasis1DQuadValues[(k % (order + 1)) * quadOrder + quadi] *
            advectionBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];

        rhoLocal += rhoK * basisk2D;
      }

      double uMeshLocal = 0.0;
      double vMeshLocal = 0.0;
      for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
        size_t indK = getKinematicIndexFromCell(celli, cellj, k);
        const double uk = uMesh(indK);
        const double vk = vMesh(indK);

        const double basisk2D =
            kinematicBasis1DQuadValues[(k % (order + 1)) * quadOrder + quadi] *
            kinematicBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];

        uMeshLocal += uk * basisk2D;
        vMeshLocal += vk * basisk2D;
      }

      for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
        const double basisiValue =
            kinematicBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                       quadi] *
            kinematicBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                       quadj];
        for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
          const double basisj2Ddx =
              kinematicBasis1DdxQuadValues[(basisj % (order + 1)) * quadOrder +
                                           quadi] *
              kinematicBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                         quadj];
          const double basisj2Ddy =
              kinematicBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                         quadi] *
              kinematicBasis1DdxQuadValues[(basisj / (order + 1)) * quadOrder +
                                           quadj];

          double dotProd = uMeshLocal * basisj2Ddx + vMeshLocal * basisj2Ddy;

          output(basisi, basisj) += legendreWeights[indMin + quadi] *
                                    legendreWeights[indMin + quadj] * rhoLocal *
                                    basisiValue * dotProd * jacDet;
        }
      }
    }
  }

  return output;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> FEMALEMethod::quadCellK(
    size_t celli, size_t cellj,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> output(
      (order + 1) * (order + 1), (order + 1) * (order + 1));

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quadi = 0; quadi < quadOrder; quadi++) {
    for (size_t quadj = 0; quadj < quadOrder; quadj++) {
      Eigen::Matrix<double, 2, 2> jacobian =
          getCellJacobian(celli, cellj, quadi, quadj, x, y);
      const double jacDet = std::abs(jacobian.determinant());

      double uMeshLocal = 0.0;
      double vMeshLocal = 0.0;
      for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
        size_t indK = getKinematicIndexFromCell(celli, cellj, k);
        const double uk = uMesh(indK);
        const double vk = vMesh(indK);

        const double basisk2D =
            kinematicBasis1DQuadValues[(k % (order + 1)) * quadOrder + quadi] *
            kinematicBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];

        uMeshLocal += uk * basisk2D;
        vMeshLocal += vk * basisk2D;
      }

      for (size_t basisi = 0; basisi < (order + 1) * (order + 1); basisi++) {
        const double basisi2D =
            advectionBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                       quadi] *
            advectionBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                       quadj];

        for (size_t basisj = 0; basisj < (order + 1) * (order + 1); basisj++) {
          const double basisj2Ddx =
              advectionBasis1DdxQuadValues[(basisj % (order + 1)) * quadOrder +
                                           quadi] *
              advectionBasis1DQuadValues[(basisj / (order + 1)) * quadOrder +
                                         quadj];
          const double basisj2Ddy =
              advectionBasis1DQuadValues[(basisj % (order + 1)) * quadOrder +
                                         quadi] *
              advectionBasis1DdxQuadValues[(basisj / (order + 1)) * quadOrder +
                                           quadj];

          double dotProd = uMeshLocal * basisj2Ddx + vMeshLocal * basisj2Ddy;

          output(basisi, basisj) += legendreWeights[indMin + quadi] *
                                    legendreWeights[indMin + quadj] * basisi2D *
                                    dotProd * jacDet;
        }
      }
    }
  }

  return output;
}

void FEMALEMethod::quadVerticalFacesK(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t facei = 1; facei < xSize; facei++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      const size_t indKinematicBottom =
          getKinematicIndexFromCell(facei, cellj, 0);
      const size_t indKinematicTop =
          getKinematicIndexFromCell(facei, cellj, order * (order + 1));
      const double dx = x(indKinematicTop) - x(indKinematicBottom);
      const double dy = y(indKinematicTop) - y(indKinematicBottom);
      const double dl = std::sqrt(dx * dx + dy * dy);
      const double nfx = dy / dl;
      const double nfy = -dx / dl;
      for (size_t quad = 0; quad < quadOrder; quad++) {
        double dXdy = 0.0;
        double dYdy = 0.0;
        double uMeshLocal = 0.0;
        double vMeshLocal = 0.0;
        for (size_t k = 0; k < order + 1; k++) {
          const size_t indK =
              getKinematicIndexFromCell(facei, cellj, k * (order + 1));
          const double xk = x(indK);
          const double yk = y(indK);
          const double uk = uMesh(indK);
          const double vk = vMesh(indK);
          const double basisk2D =
              advectionBasis1DQuadValues[k * quadOrder + quad];
          const double basisk2Ddy =
              advectionBasis1DdxQuadValues[k * quadOrder + quad];

          dXdy += xk * basisk2Ddy;
          dYdy += yk * basisk2Ddy;
          uMeshLocal += uk * basisk2D;
          vMeshLocal += vk * basisk2D;
        }
        const double localScale = std::sqrt(dXdy * dXdy + dYdy * dYdy);
        const double dotProd = uMeshLocal * nfx + vMeshLocal * nfy;

        for (size_t basisi = 0; basisi < order + 1; basisi++) {
          const double basisi2D =
              advectionBasis1DQuadValues[basisi * quadOrder + quad];
          const double flux = dotProd * basisi2D;
          size_t indI = 0;
          if (dotProd >= 0.0) {
            indI = getAdvectionIndexFromCell(facei - 1, cellj,
                                             basisi * (order + 1) + order);
          } else {
            indI =
                getAdvectionIndexFromCell(facei, cellj, basisi * (order + 1));
          }
          for (size_t basisj = 0; basisj < order + 1; basisj++) {
            const double basisj2D =
                advectionBasis1DQuadValues[basisj * quadOrder + quad];
            const double quadValue =
                legendreWeights[indMin + quad] * flux * basisj2D * localScale;
            const size_t indJleft = getAdvectionIndexFromCell(
                facei - 1, cellj, basisj * (order + 1) + order);
            const size_t indJright =
                getAdvectionIndexFromCell(facei, cellj, basisj * (order + 1));
            K.coeffRef(indI, indJleft) -= quadValue;
            K.coeffRef(indI, indJright) += quadValue;
          }
        }
      }
    }
  }
}

void FEMALEMethod::quadHorizontalFacesK(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t facej = 1; facej < ySize; facej++) {
      const size_t indKinematicLeft =
          getKinematicIndexFromCell(celli, facej, 0);
      const size_t indKinematicRight =
          getKinematicIndexFromCell(celli, facej, order);
      const double dx = x(indKinematicRight) - x(indKinematicLeft);
      const double dy = y(indKinematicRight) - y(indKinematicLeft);
      const double dl = std::sqrt(dx * dx + dy * dy);
      const double nfx = -dy / dl;
      const double nfy = dx / dl;
      for (size_t quad = 0; quad < quadOrder; quad++) {
        double dXdx = 0.0;
        double dYdx = 0.0;
        double uMeshLocal = 0.0;
        double vMeshLocal = 0.0;
        for (size_t k = 0; k < order + 1; k++) {
          const size_t indK = getKinematicIndexFromCell(celli, facej, k);
          const double xk = x(indK);
          const double yk = y(indK);
          const double uk = uMesh(indK);
          const double vk = vMesh(indK);
          const double basisk2D =
              advectionBasis1DQuadValues[k * quadOrder + quad];
          const double basisk2Ddx =
              advectionBasis1DdxQuadValues[k * quadOrder + quad];

          dXdx += xk * basisk2Ddx;
          dYdx += yk * basisk2Ddx;
          uMeshLocal += uk * basisk2D;
          vMeshLocal += vk * basisk2D;
        }
        const double localScale = std::sqrt(dXdx * dXdx + dYdx * dYdx);
        const double dotProd = uMeshLocal * nfx + vMeshLocal * nfy;

        for (size_t basisi = 0; basisi < order + 1; basisi++) {
          const double basisi2D =
              advectionBasis1DQuadValues[basisi * quadOrder + quad];
          const double flux = dotProd * basisi2D;
          size_t indI = 0;
          if (dotProd >= 0.0) {
            indI = getAdvectionIndexFromCell(celli, facej - 1,
                                             order * (order + 1) + basisi);
          } else {
            indI = getAdvectionIndexFromCell(celli, facej, basisi);
          }
          for (size_t basisj = 0; basisj < order + 1; basisj++) {
            const double basisj2D =
                advectionBasis1DQuadValues[basisj * quadOrder + quad];
            const double quadValue =
                legendreWeights[indMin + quad] * flux * basisj2D * localScale;
            const size_t indJdown = getAdvectionIndexFromCell(
                celli, facej - 1, order * (order + 1) + basisj);
            const size_t indJup =
                getAdvectionIndexFromCell(celli, facej, basisj);
            K.coeffRef(indI, indJdown) -= quadValue;
            K.coeffRef(indI, indJup) += quadValue;
          }
        }
      }
    }
  }
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

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t i = 0; i < quadOrder; i++) {
    for (size_t j = 0; j < quadOrder; j++) {
      Eigen::Matrix<double, 2, 2> jacobian =
          getCellJacobian(celli, cellj, i, j, x, y);
      Eigen::JacobiSVD<Eigen::Matrix<double, 2, 2>> svd(jacobian);
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
              kinematicBasis1DdxQuadValues[(basisKinematic % (order + 1)) *
                                               quadOrder +
                                           i] *
              kinematicBasis1DQuadValues[(basisKinematic / (order + 1)) *
                                             quadOrder +
                                         j];
          double basisKinematic2Ddy =
              kinematicBasis1DQuadValues[(basisKinematic % (order + 1)) *
                                             quadOrder +
                                         i] *
              kinematicBasis1DdxQuadValues[(basisKinematic / (order + 1)) *
                                               quadOrder +
                                           j];
          Eigen::Matrix<double, 2, 2> gradBasisx{
              {basisKinematic2Ddx, basisKinematic2Ddy}, {0.0, 0.0}};
          Eigen::Matrix<double, 2, 2> gradBasisy{
              {0.0, 0.0}, {basisKinematic2Ddx, basisKinematic2Ddy}};
          double thermobasis =
              thermoBasis1DQuadValues[(basisThermo % order) * quadOrder + i] *
              thermoBasis1DQuadValues[(basisThermo / order) * quadOrder + j];

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
size_t FEMALEMethod::getAdvectionIndexFromCell(const size_t celli,
                                               const size_t cellj,
                                               const size_t k) const {
  assert(k < (order + 1) * (order + 1));
  return (celli * (order + 1) * (order + 1) * ySize +
          cellj * (order + 1) * (order + 1) + k);
}
size_t FEMALEMethod::getQuadIndexFromCell(const size_t celli,
                                          const size_t cellj,
                                          const size_t k) const {
  assert(k < quadOrder * quadOrder);
  return celli * quadOrder * quadOrder * ySize + cellj * quadOrder * quadOrder +
         k;
}

double FEMALEMethod::getMinRhoFromCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double rhoMin = std::numeric_limits<double>::max();
  for (size_t k = 0; k < quadOrder * quadOrder; k++) {
    const size_t indRho = getQuadIndexFromCell(celli, cellj, k);
    const double rho = rhoQuad(indRho);
    if (rho < rhoMin) {
      rhoMin = rho;
    }
  }
  return rhoMin;
}

double FEMALEMethod::getMinRhoFromCellAndNeighbours(size_t celli,
                                                    size_t cellj) {
  double rhoMin = getMinRhoFromCell(celli, cellj);
  if (celli > 0) {
    const double cellRhoMin = getMinRhoFromCell(celli - 1, cellj);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  if (celli + 1 < xSize) {
    const double cellRhoMin = getMinRhoFromCell(celli + 1, cellj);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  if (cellj > 0) {
    const double cellRhoMin = getMinRhoFromCell(celli, cellj - 1);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellRhoMin = getMinRhoFromCell(celli, cellj + 1);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  return rhoMin;
}

double FEMALEMethod::getMaxRhoFromCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double rhoMax = std::numeric_limits<double>::min();
  for (size_t k = 0; k < quadOrder * quadOrder; k++) {
    const size_t indRho = getQuadIndexFromCell(celli, cellj, k);
    const double rho = rhoQuad(indRho);
    if (rho > rhoMax) {
      rhoMax = rho;
    }
  }
  return rhoMax;
}

double FEMALEMethod::getMaxRhoFromCellAndNeighbours(size_t celli,
                                                    size_t cellj) {
  double rhoMax = getMaxRhoFromCell(celli, cellj);
  if (celli > 0) {
    const double cellRhoMax = getMaxRhoFromCell(celli - 1, cellj);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  if (celli + 1 < xSize) {
    const double cellRhoMax = getMaxRhoFromCell(celli + 1, cellj);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  if (cellj > 0) {
    const double cellRhoMax = getMaxRhoFromCell(celli, cellj - 1);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellRhoMax = getMaxRhoFromCell(celli, cellj + 1);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  return rhoMax;
}

double FEMALEMethod::getMinEFromCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double eMin = std::numeric_limits<double>::max();
  for (size_t k = 0; k < order * order; k++) {
    const size_t indK = getThermodynamicIndexFromCell(celli, cellj, k);
    const double eK = e(indK);
    if (eK < eMin) {
      eMin = eK;
    }
  }
  return eMin;
}

double FEMALEMethod::getMinEFromCellAndNeighbours(size_t celli, size_t cellj) {
  double eMin = getMinEFromCell(celli, cellj);
  if (celli > 0) {
    const double cellEMin = getMinEFromCell(celli - 1, cellj);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  if (celli + 1 < xSize) {
    const double cellEMin = getMinEFromCell(celli + 1, cellj);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  if (cellj > 0) {
    const double cellEMin = getMinEFromCell(celli, cellj - 1);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellEMin = getMinEFromCell(celli, cellj + 1);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  return eMin;
}

double FEMALEMethod::getMaxEFromCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double eMax = std::numeric_limits<double>::min();
  for (size_t k = 0; k < order * order; k++) {
    const size_t indK = getThermodynamicIndexFromCell(celli, cellj, k);
    const double eK = e(indK);
    if (eK > eMax) {
      eMax = eK;
    }
  }
  return eMax;
}

double FEMALEMethod::getMaxEFromCellAndNeighbours(size_t celli, size_t cellj) {
  double eMax = getMaxEFromCell(celli, cellj);
  if (celli > 0) {
    const double cellEMax = getMaxEFromCell(celli - 1, cellj);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  if (celli + 1 < xSize) {
    const double cellEMax = getMaxEFromCell(celli + 1, cellj);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  if (cellj > 0) {
    const double cellEMax = getMaxEFromCell(celli, cellj - 1);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellEMax = getMaxEFromCell(celli, cellj + 1);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  return eMax;
}

double FEMALEMethod::getMinRhoFromRemapCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double rhoMin = std::numeric_limits<double>::max();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
    const double rhoLocal = rhoRemap(indK);
    if (rhoLocal < rhoMin) {
      rhoMin = rhoLocal;
    }
  }
  return rhoMin;
}

double FEMALEMethod::getMinRhoFromRemapCellAndNeighbours(size_t celli,
                                                         size_t cellj) {
  double rhoMin = getMinRhoFromRemapCell(celli, cellj);
  if (celli > 0) {
    const double cellRhoMin = getMinRhoFromRemapCell(celli - 1, cellj);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  if (celli + 1 < xSize) {
    const double cellRhoMin = getMinRhoFromRemapCell(celli + 1, cellj);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  if (cellj > 0) {
    const double cellRhoMin = getMinRhoFromRemapCell(celli, cellj - 1);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellRhoMin = getMinRhoFromRemapCell(celli, cellj + 1);
    if (cellRhoMin < rhoMin) {
      rhoMin = cellRhoMin;
    }
  }
  return rhoMin;
}

double FEMALEMethod::getMaxRhoFromRemapCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double rhoMax = std::numeric_limits<double>::min();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
    const double rhoLocal = rhoRemap(indK);
    if (rhoLocal > rhoMax) {
      rhoMax = rhoLocal;
    }
  }
  return rhoMax;
}

double FEMALEMethod::getMaxRhoFromRemapCellAndNeighbours(size_t celli,
                                                         size_t cellj) {
  double rhoMax = getMaxRhoFromRemapCell(celli, cellj);
  if (celli > 0) {
    const double cellRhoMax = getMaxRhoFromRemapCell(celli - 1, cellj);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  if (celli + 1 < xSize) {
    const double cellRhoMax = getMaxRhoFromRemapCell(celli + 1, cellj);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  if (cellj > 0) {
    const double cellRhoMax = getMaxRhoFromRemapCell(celli, cellj - 1);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellRhoMax = getMaxRhoFromRemapCell(celli, cellj + 1);
    if (cellRhoMax > rhoMax) {
      rhoMax = cellRhoMax;
    }
  }
  return rhoMax;
}

double FEMALEMethod::getMinEFromRemapCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double eMin = std::numeric_limits<double>::max();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
    const double rhoEK = rhoERemap(indK);
    const double rhoK = rhoRemap(indK);
    const double eK = rhoEK / rhoK;
    if (eK < eMin) {
      eMin = eK;
    }
  }
  return eMin;
}

double FEMALEMethod::getMinEFromRemapCellAndNeighbours(size_t celli,
                                                       size_t cellj) {
  double eMin = getMinEFromRemapCell(celli, cellj);
  if (celli > 0) {
    const double cellEMin = getMinEFromRemapCell(celli - 1, cellj);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  if (celli + 1 < xSize) {
    const double cellEMin = getMinEFromRemapCell(celli + 1, cellj);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  if (cellj > 0) {
    const double cellEMin = getMinEFromRemapCell(celli, cellj - 1);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellEMin = getMinEFromRemapCell(celli, cellj + 1);
    if (cellEMin < eMin) {
      eMin = cellEMin;
    }
  }
  return eMin;
}

double FEMALEMethod::getMaxEFromRemapCell(size_t celli, size_t cellj) {
  assert(celli < xSize);
  assert(cellj < ySize);
  double eMax = std::numeric_limits<double>::min();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
    const double rhoEK = rhoERemap(indK);
    const double rhoK = rhoRemap(indK);
    const double eK = rhoEK / rhoK;
    if (eK > eMax) {
      eMax = eK;
    }
  }
  return eMax;
}
double FEMALEMethod::getMaxEFromRemapCellAndNeighbours(size_t celli,
                                                       size_t cellj) {
  double eMax = getMaxEFromRemapCell(celli, cellj);
  if (celli > 0) {
    const double cellEMax = getMaxEFromRemapCell(celli - 1, cellj);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  if (celli + 1 < xSize) {
    const double cellEMax = getMaxEFromRemapCell(celli + 1, cellj);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  if (cellj > 0) {
    const double cellEMax = getMaxEFromRemapCell(celli, cellj - 1);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  if (cellj + 1 < ySize) {
    const double cellEMax = getMaxEFromRemapCell(celli, cellj + 1);
    if (cellEMax > eMax) {
      eMax = cellEMax;
    }
  }
  return eMax;
}

void FEMALEMethod::initInitializers() {
  // rhoInitializer =
  //     [xmin = problem.xmin, xmax = problem.xmax, ymin = problem.ymin,
  //      ymax = problem.ymax, xSize = this->xSize, ySize = this->ySize,
  //      rhoInitializer = problem.rhoInitializer](double x, double y) {
  //       double dx = (xmax - xmin) / xSize;
  //       double dy = (ymax - ymin) / ySize;
  //       double xij = xmin + (std::floor((x - xmin) / dx) + 0.5) * dx;
  //       double yij = ymin + (std::floor((y - ymin) / dy) + 0.5) * dy;
  //       return rhoInitializer(xij, yij);
  //     };
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
  initQuadBasisValues();
  initOutputBasisValues();
}

void FEMALEMethod::initQuadBasisValues() {
  kinematicBasis1DQuadValues.reserve((order + 1) * quadOrder);
  kinematicBasis1DdxQuadValues.reserve((order + 1) * quadOrder);

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t k = 0; k < order + 1; k++) {
    for (size_t i = indMin; i < indMin + quadOrder; i++) {
      double xi = legendreAbscissas[i];
      kinematicBasis1DQuadValues.push_back(lobattoBasis1D(xi, order, k));
      kinematicBasis1DdxQuadValues.push_back(lobattoBasis1Ddx(xi, order, k));
    }
  }

  thermoBasis1DQuadValues.reserve(order * quadOrder);
  for (size_t k = 0; k < order; k++) {
    for (size_t i = indMin; i < indMin + quadOrder; i++) {
      double xi = legendreAbscissas[i];
      thermoBasis1DQuadValues.push_back(legendreBasis1D(xi, order - 1, k));
    }
  }

  advectionBasis1DQuadValues.reserve((order + 1) * quadOrder);
  advectionBasis1DdxQuadValues.reserve((order + 1) * quadOrder);
  for (size_t k = 0; k < order + 1; k++) {
    for (size_t i = indMin; i < indMin + quadOrder; i++) {
      double xi = legendreAbscissas[i];
      advectionBasis1DQuadValues.push_back(bernsteinBasis1D(xi, order, k));
      advectionBasis1DdxQuadValues.push_back(bernsteinBasis1Ddx(xi, order, k));
    }
  }
}

void FEMALEMethod::initOutputBasisValues() {
  output1DKinematicValues.reserve((order + 1) * order);
  output1DdxKinematicValues.reserve((order + 1) * order);
  output1DThermoValues.reserve(order * order);
  output1DRhoValues.reserve(quadOrder * order);

  const double dx = 1.0 / order;
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
      output1DThermoValues.push_back(legendreBasis1D(xi, order - 1, k));
    }
  }

  for (size_t k = 0; k < quadOrder; k++) {
    for (size_t i = 0; i < order; i++) {
      double xi = 0.5 * dx + dx * i;
      output1DRhoValues.push_back(legendreBasis1D(xi, quadOrder - 1, k));
    }
  }
}

void FEMALEMethod::initRhoValues() {
  const double celldx = (problem.xmax - problem.xmin) / xSize;
  const double celldy = (problem.ymax - problem.ymin) / ySize;
  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t k = 0; k < quadOrder * quadOrder; k++) {
        double xlocal = legendreAbscissas[indMin + k % quadOrder];
        double ylocal = legendreAbscissas[indMin + k / quadOrder];
        double xij = problem.xmin + celldx * (celli + xlocal);
        double yij = problem.ymin + celldy * (cellj + ylocal);

        const size_t indRho = getQuadIndexFromCell(celli, cellj, k);
        rhoQuad(indRho) = problem.rhoInitializer(xij, yij);
      }
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
        double xij = problem.xmin + celldx * (celli + localdx);
        double yij = problem.ymin + celldy * (cellj + localdy);
        auto eos = eosInitializer(xij, yij);
        auto p = problem.pInitializer(xij, yij);
        auto rho = problem.rhoInitializer(xij, yij);
        e(getThermodynamicIndexFromCell(celli, cellj, k)) = eos->gete(rho, p);
      }
    }
  }
}

void FEMALEMethod::initKinematicMassMatrix() { calcKinematicMassMatrix(); }

void FEMALEMethod::initThermodynamicInverseMassMatrix() {
  Mt_inv.reserve(Eigen::VectorXi::Constant(Nt, order * order));
  calcThermoMassMatrix();
  Mt_inv.makeCompressed();
}

void FEMALEMethod::initForceMatrices() {
  Fx.reserve(Eigen::VectorXi::Constant(Nt, (order + 1) * (order + 1)));
  Fy.reserve(Eigen::VectorXi::Constant(Nt, (order + 1) * (order + 1)));
}

void FEMALEMethod::initVectorMatrices() {
  vectorLaplacianX.reserve(
      Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));
  vectorLaplacianY.reserve(
      Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));
  vectorMass.reserve(
      Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));

  const Eigen::Matrix<double, 2, 2> cellVectorMass = quadCellVectorMass();
  const Eigen::Matrix<double, 2, 2> cellVectorLaplacian =
      quadCellVectorLaplacian();

  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
          vectorMass.coeffRef(indI, indJ) += cellVectorMass(i, j);
          vectorLaplacianX.coeffRef(indI, indJ) += cellVectorLaplacian(i, j);
          vectorLaplacianY.coeffRef(indI, indJ) += cellVectorLaplacian(i, j);
        }
      }
    }
  }
  resolveBoundaryVector();
  vectorMass.makeCompressed();
  vectorLaplacianX.makeCompressed();
  vectorLaplacianY.makeCompressed();
}

void FEMALEMethod::initRemapMatrices() {
  Mv.reserve(Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));
  M.reserve(Eigen::VectorXi::Constant(Na, (2 * order + 1) * (2 * order + 1)));
  Kv.reserve(Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));
  MInv.reserve(
      Eigen::VectorXi::Constant(Na, (2 * order + 1) * (2 * order + 1)));
}

void FEMALEMethod::initSolvers() {
  kinematicSolverx.compute(Mkx);
  assert(kinematicSolverx.info() == Eigen::Success);
  kinematicSolvery.compute(Mky);
  assert(kinematicSolvery.info() == Eigen::Success);
}

void FEMALEMethod::calcKinematicMassMatrix() {
  Mkx.setZero();
  Mkx.reserve(Eigen::VectorXi::Constant(
      Nk, (2 * order + 1) *
              (2 * order + 1)));  // reserve maximum possible space for each row
  Mky.setZero();
  Mky.reserve(Eigen::VectorXi::Constant(
      Nk, (2 * order + 1) *
              (2 * order + 1)));  // reserve maximum possible space for each row
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> quadResult =
          quadKinematicCellMass(celli, cellj);

      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          const size_t indJ = getKinematicIndexFromCell(celli, cellj, j);

          Mkx.coeffRef(indI, indJ) += quadResult(i, j);
          Mky.coeffRef(indI, indJ) += quadResult(i, j);
        }
      }
    }
  }
  resolveBoundaryMass();
  Mkx.makeCompressed();
  Mky.makeCompressed();
}

void FEMALEMethod::calcThermoMassMatrix() {
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MtCell =
          quadThermoCellMass(celli, cellj);

      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MtCellInv =
          MtCell.inverse();

      for (size_t i = 0; i < order * order; i++) {
        const size_t indI = getThermodynamicIndexFromCell(celli, cellj, i);
        for (size_t j = 0; j < order * order; j++) {
          const size_t indJ = getThermodynamicIndexFromCell(celli, cellj, j);
          Mt_inv.coeffRef(indI, indJ) = MtCellInv(i, j);
        }
      }
    }
  }
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

  const size_t indMin = getLegendreStartIndex(quadOrder);
  const double celldx = (problem.xmax - problem.xmin) / xSize;
  const double celldy = (problem.ymax - problem.ymin) / ySize;
  const double xlocal = legendreAbscissas[indMin + i];
  const double ylocal = legendreAbscissas[indMin + j];
  const double xij = problem.xmin + celldx * (celli + xlocal);
  const double yij = problem.ymin + celldy * (cellj + ylocal);

  double eLocal = 0.0;
  for (size_t k = 0; k < order * order; k++) {
    size_t indK = getThermodynamicIndexFromCell(celli, cellj, k);
    eLocal += e[indK] * thermoBasis1DQuadValues[(k % order) * quadOrder + i] *
              thermoBasis1DQuadValues[(k / order) * quadOrder + j];
  }

  Eigen::Matrix<double, 2, 2> jacobianInitial =
      getCellJacobian(celli, cellj, i, j, xInitial, yInitial);

  // double rhoInitial = rhoInitializer(xij, yij);
  const size_t indRho = getQuadIndexFromCell(celli, cellj, i * quadOrder + j);
  double rhoInitial = rhoQuad(indRho);
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
        kinematicBasis1DdxQuadValues[(k % (order + 1)) * quadOrder + i] *
        kinematicBasis1DQuadValues[(k / (order + 1)) * quadOrder + j];
    double basis2Ddy =
        kinematicBasis1DQuadValues[(k % (order + 1)) * quadOrder + i] *
        kinematicBasis1DdxQuadValues[(k / (order + 1)) * quadOrder + j];
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

void FEMALEMethod::calcRemapMatrices(
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  calcMv(x, y);
  calcM(x, y);
  calcKv(x, y);
  calcK(x, y);
  calcAuxMat();
  calcRemapTau();
}

void FEMALEMethod::calcMv(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  Mv.setZero();
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellMv =
          quadCellMv(celli, cellj, x, y);

      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          size_t indJ = getKinematicIndexFromCell(celli, cellj, j);

          Mv.coeffRef(indI, indJ) += cellMv(i, j);
        }
      }
    }
  }

  if (!Mv.isCompressed()) {
    Mv.makeCompressed();
  }
  MvSolver.compute(Mv);
}

void FEMALEMethod::calcM(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellM =
          quadCellM(celli, cellj, x, y);
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellMInv =
          cellM.inverse();

      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
          size_t indJ = getAdvectionIndexFromCell(celli, cellj, j);

          M.coeffRef(indI, indJ) = cellM(i, j);
          MInv.coeffRef(indI, indJ) = cellMInv(i, j);
        }
      }
    }
  }

  if (!M.isCompressed()) {
    M.makeCompressed();
  }
  if (!MInv.isCompressed()) {
    MInv.makeCompressed();
  }
}
void FEMALEMethod::calcKv(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                          const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  Kv.setZero();
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellKv =
          quadCellKv(celli, cellj, x, y);

      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          size_t indJ = getKinematicIndexFromCell(celli, cellj, j);

          Kv.coeffRef(indI, indJ) += cellKv(i, j);
        }
      }
    }
  }

  if (!Kv.isCompressed()) {
    Kv.makeCompressed();
  }
}
void FEMALEMethod::calcK(const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
                         const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) {
  K.setZero();
  K.reserve(Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> cellK =
          quadCellK(celli, cellj, x, y);

      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
          size_t indJ = getAdvectionIndexFromCell(celli, cellj, j);

          K.coeffRef(indI, indJ) = cellK(i, j);
        }
      }
    }
  }

  quadVerticalFacesK(x, y);
  quadHorizontalFacesK(x, y);

  K.makeCompressed();
}

void FEMALEMethod::calcAuxMat() {
  // L and Mstar;
  L = -M;
  for (size_t i = 0; i < Na; i++) {
    L.coeffRef(i, i) = 0.0;
    for (size_t j = 0; j < i; j++) {
      L.coeffRef(i, i) -= L.coeff(i, j);
    }
    for (size_t j = i + 1; j < Na; j++) {
      L.coeffRef(i, i) -= L.coeff(i, j);
    }
  }
  Mlumped = M + L;
  for (size_t i = 0; i < Na; i++) {
    const double Mlumpedii = Mlumped.coeff(i, i);
    assert(Mlumpedii != 0.0);
    MlumpedInv.coeffRef(i, i) = 1.0 / Mlumpedii;
  }
  if (!MlumpedInv.isCompressed()) {
    MlumpedInv.makeCompressed();
  }

  // D and Kstar;
  D.setZero();
  D.reserve(Eigen::VectorXi::Constant(Nk, (2 * order + 1) * (2 * order + 1)));
  for (size_t i = 0; i < Na; i++) {
    for (size_t j = 0; j < i; j++) {
      double Dij = std::max(-K.coeff(i, j), -K.coeff(j, i));
      Dij = std::max(0.0, Dij);
      if (Dij == 0.0) {
        continue;
      }
      D.coeffRef(i, j) = Dij;
      D.coeffRef(i, i) -= Dij;
    }
    for (size_t j = i + 1; j < Na; j++) {
      double Dij = std::max(-K.coeff(i, j), -K.coeff(j, i));
      Dij = std::max(0.0, Dij);
      if (Dij == 0.0) {
        continue;
      }
      D.coeffRef(i, j) = Dij;
      D.coeffRef(i, i) -= Dij;
    }
  }
  D.makeCompressed();
  Kupwinded = K + D;
}

void FEMALEMethod::optimizeMesh() {
  xOptimal = x;
  yOptimal = y;
  Eigen::Matrix<double, Eigen::Dynamic, 1> xTmp(Nk);
  Eigen::Matrix<double, Eigen::Dynamic, 1> yTmp(Nk);

  // DECOMPOSITION A = D + L + U
  // TEST with 0.5 * eps to account for multiple dimensions
  Eigen::SparseMatrix<double> Ax = eps * vectorLaplacianX + vectorMass;
  Eigen::SparseMatrix<double> Ay = eps * vectorLaplacianY + vectorMass;
  Ax.makeCompressed();
  Ay.makeCompressed();

  Eigen::SparseMatrix<double> Dx_inv(Nk, Nk);
  Eigen::SparseMatrix<double> Dy_inv(Nk, Nk);
  Dx_inv.reserve(Eigen::VectorXi::Constant(Nk, 1));
  Dy_inv.reserve(Eigen::VectorXi::Constant(Nk, 1));
  for (size_t i = 0; i < Nk; i++) {
    double Axii = Ax.coeff(i, i);
    double Ayii = Ay.coeff(i, i);
    Dx_inv.insert(i, i) = 1.0 / Axii;
    Dy_inv.insert(i, i) = 1.0 / Ayii;
  }

  Eigen::SparseMatrix<double> LUx = Ax;
  Eigen::SparseMatrix<double> LUy = Ay;
  for (size_t i = 0; i < Nk; i++) {
    LUx.coeffRef(i, i) = 0.0;
    LUy.coeffRef(i, i) = 0.0;
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> bx = vectorMass * x;
  Eigen::Matrix<double, Eigen::Dynamic, 1> by = vectorMass * y;

  constexpr size_t nSteps = 10;
  for (size_t iter = 0; iter < nSteps; iter++) {
    xTmp.swap(xOptimal);
    yTmp.swap(yOptimal);
    xOptimal = Dx_inv * (bx - LUx * xTmp);
    yOptimal = Dy_inv * (by - LUy * yTmp);
  }

  uMesh = xOptimal - x;
  vMesh = yOptimal - y;
}

void FEMALEMethod::transitionToRemap() {
  preTransitionRhoToRemap();
  transitionRhoToRemap();
  transitionEToRemap();
}

void FEMALEMethod::preTransitionRhoToRemap() {
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
        const size_t quadi = quad % quadOrder;
        const size_t quadj = quad / quadOrder;

        Eigen::Matrix<double, 2, 2> jacobian =
            getCellJacobian(celli, cellj, quadi, quadj, x, y);
        const double jacDet = std::abs(jacobian.determinant());
        jacobian =
            getCellJacobian(celli, cellj, quadi, quadj, xInitial, yInitial);
        const double jacDetInitial = std::abs(jacobian.determinant());

        const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
        rhoQuad(indRho) *= jacDetInitial / jacDet;
      }
    }
  }
}

void FEMALEMethod::transitionRhoToRemap() {
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      preRhoToRemap(celli, cellj);
      FCTProjection();
      postRhoToRemap(celli, cellj);
    }
  }
}

void FEMALEMethod::preRhoToRemap(size_t celli, size_t cellj) {
  const size_t nProj = (order + 1) * (order + 1);

  wQuad.resize(quadOrder * quadOrder);
  detJacQuad.resize(quadOrder * quadOrder);
  basisQuadValues.resize(nProj, quadOrder * quadOrder);
  xQuad.resize(quadOrder * quadOrder);
  xProj.resize(nProj);
  yQuad.resize(quadOrder * quadOrder);
  ymin.resize(nProj);
  ymax.resize(nProj);

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const size_t quadi = quad % quadOrder;
    const size_t quadj = quad / quadOrder;

    wQuad(quad) =
        legendreWeights[indMin + quadi] * legendreWeights[indMin + quadj];

    Eigen::Matrix<double, 2, 2> jacobian =
        getCellJacobian(celli, cellj, quadi, quadj, x, y);
    detJacQuad(quad) = std::abs(jacobian.determinant());

    const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
    yQuad(quad) = rhoQuad(indRho);
  }
  xQuad = Eigen::VectorXd::Constant(quadOrder * quadOrder, 1.0);
  xProj = Eigen::VectorXd::Constant(nProj, 1.0);
  for (size_t basisi = 0; basisi < nProj; basisi++) {
    for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
      const size_t quadi = quad % quadOrder;
      const size_t quadj = quad / quadOrder;
      basisQuadValues(basisi, quad) =
          advectionBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                     quadi] *
          advectionBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                     quadj];
    }
  }
  double rhoMin = getMinRhoFromCellAndNeighbours(celli, cellj);
  double rhoMax = getMaxRhoFromCellAndNeighbours(celli, cellj);
  ymin = Eigen::VectorXd::Constant(nProj, rhoMin);
  ymax = Eigen::VectorXd::Constant(nProj, rhoMax);
}

void FEMALEMethod::postRhoToRemap(size_t celli, size_t cellj) {
  const size_t nProj = (order + 1) * (order + 1);
  assert(xyOut.size() == nProj);
  for (size_t k = 0; k < nProj; k++) {
    const size_t indRho = getAdvectionIndexFromCell(celli, cellj, k);
    rhoRemap(indRho) = xyOut(k);
  }
}

void FEMALEMethod::transitionEToRemap() {
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      preEToRemap(celli, cellj);
      FCTProjection();
      postEToRemap(celli, cellj);
    }
  }
}

void FEMALEMethod::preEToRemap(size_t celli, size_t cellj) {
  const size_t nProj = (order + 1) * (order + 1);

  wQuad.resize(quadOrder * quadOrder);
  detJacQuad.resize(quadOrder * quadOrder);
  basisQuadValues.resize(nProj, quadOrder * quadOrder);
  xQuad.resize(quadOrder * quadOrder);
  xProj.resize(nProj);
  yQuad.resize(quadOrder * quadOrder);
  ymin.resize(nProj);
  ymax.resize(nProj);

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const size_t quadi = quad % quadOrder;
    const size_t quadj = quad / quadOrder;

    wQuad(quad) =
        legendreWeights[indMin + quadi] * legendreWeights[indMin + quadj];

    Eigen::Matrix<double, 2, 2> jacobian =
        getCellJacobian(celli, cellj, quadi, quadj, x, y);
    detJacQuad(quad) = std::abs(jacobian.determinant());

    const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
    xQuad(quad) = rhoQuad(indRho);

    double eQuad = 0.0;
    for (size_t k = 0; k < order * order; k++) {
      const size_t indK = getThermodynamicIndexFromCell(celli, cellj, k);
      const double ek = e(indK);
      const double basisk2D =
          thermoBasis1DQuadValues[(k % order) * quadOrder + quadi] *
          thermoBasis1DQuadValues[(k / order) * quadOrder + quadj];

      eQuad += ek * basisk2D;
    }
    yQuad(quad) = eQuad;
  }
  for (size_t basisi = 0; basisi < nProj; basisi++) {
    const size_t indRho = getAdvectionIndexFromCell(celli, cellj, basisi);
    xProj(basisi) = rhoRemap(indRho);
  }
  for (size_t basisi = 0; basisi < nProj; basisi++) {
    for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
      const size_t quadi = quad % quadOrder;
      const size_t quadj = quad / quadOrder;
      basisQuadValues(basisi, quad) =
          advectionBasis1DQuadValues[(basisi % (order + 1)) * quadOrder +
                                     quadi] *
          advectionBasis1DQuadValues[(basisi / (order + 1)) * quadOrder +
                                     quadj];
    }
  }

  double eMin = getMinEFromCellAndNeighbours(celli, cellj);
  double eMax = getMaxEFromCellAndNeighbours(celli, cellj);
  ymin = Eigen::VectorXd::Constant(nProj, eMin);
  ymax = Eigen::VectorXd::Constant(nProj, eMax);
}

void FEMALEMethod::postEToRemap(size_t celli, size_t cellj) {
  const size_t nProj = (order + 1) * (order + 1);
  assert(xyOut.size() == nProj);
  for (size_t k = 0; k < nProj; k++) {
    const size_t indRhoE = getAdvectionIndexFromCell(celli, cellj, k);
    rhoERemap(indRhoE) = xyOut(k);
  }
}

void FEMALEMethod::transitionToLagrange() {
  transitionRhoToLagrange();
  transitionEToLagrange();
  calcKinematicMassMatrix();
  calcThermoMassMatrix();
  xInitial = x;
  yInitial = y;
}

void FEMALEMethod::transitionRhoToLagrange() {
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
        const size_t quadi = quad % quadOrder;
        const size_t quadj = quad / quadOrder;

        double rhoLocal = 0.0;
        for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
          const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
          const double rhoK = rhoRemap(indK);
          const double basisk2D =
              advectionBasis1DQuadValues[(k % (order + 1)) * quadOrder +
                                         quadi] *
              advectionBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];

          rhoLocal += rhoK * basisk2D;
        }

        const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
        rhoQuad(indRho) = rhoLocal;
      }
    }
  }
}

void FEMALEMethod::transitionEToLagrange() {}

void FEMALEMethod::preEToLagrange(size_t celli, size_t cellj) {
  const size_t nProj = order * order;

  wQuad.resize(quadOrder * quadOrder);
  detJacQuad.resize(quadOrder * quadOrder);
  basisQuadValues.resize(nProj, quadOrder * quadOrder);
  xQuad.resize(quadOrder * quadOrder);
  xProj.resize(nProj);
  yQuad.resize(quadOrder * quadOrder);
  ymin.resize(nProj);
  ymax.resize(nProj);

  const size_t indMin = getLegendreStartIndex(quadOrder);
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const size_t quadi = quad % quadOrder;
    const size_t quadj = quad / quadOrder;

    wQuad(quad) =
        legendreWeights[indMin + quadi] * legendreWeights[indMin + quadj];

    Eigen::Matrix<double, 2, 2> jacobian =
        getCellJacobian(celli, cellj, quadi, quadj, x, y);
    const size_t indRho = getQuadIndexFromCell(celli, cellj, quad);
    const double rhoLocal = rhoQuad(indRho);
    detJacQuad(quad) = std::abs(jacobian.determinant()) * rhoLocal;

    double eQuad = 0.0;
    for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
      const size_t indK = getAdvectionIndexFromCell(celli, cellj, k);
      const double eK = rhoERemap(indK);
      const double basisk2D =
          advectionBasis1DQuadValues[(k % (order + 1)) * quadOrder + quadi] *
          advectionBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];
      eQuad += eK * basisk2D;
    }
    yQuad(quad) = eQuad / rhoLocal;
  }
  xQuad = Eigen::VectorXd::Constant(quadOrder * quadOrder, 1.0);
  xProj = Eigen::VectorXd::Constant(nProj, 1.0);
  for (size_t basisi = 0; basisi < nProj; basisi++) {
    for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
      const size_t quadi = quad % quadOrder;
      const size_t quadj = quad / quadOrder;
      basisQuadValues(basisi, quad) =
          thermoBasis1DQuadValues[(basisi % order) * quadOrder + quadi] *
          thermoBasis1DQuadValues[(basisi / order) * quadOrder + quadj];
    }
  }

  double eMin = getMinEFromRemapCellAndNeighbours(celli, cellj);
  double eMax = getMaxEFromRemapCellAndNeighbours(celli, cellj);
  ymin = Eigen::VectorXd::Constant(nProj, eMin);
  ymax = Eigen::VectorXd::Constant(nProj, eMax);
}

void postEToLagrange(size_t celli, size_t cellj);

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

void FEMALEMethod::calcRemapTau() {
  for (size_t i = 0; i < Na; i++) {
    const double Mii = Mlumped.coeff(i, i);
    const double Kii = Kupwinded.coeff(i, i);
    assert(Mii != 0.0);
    assert(Kii != 0.0);
    const double dtTmp = remapCFL * Mii / Kii;
    if (dtTmp < remapTau) {
      remapTau = dtTmp;
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

void FEMALEMethod::resolveBoundaryVector() {
  resolveLeftBoundaryVector();
  resolveTopBoundaryVector();
  resolveRightBoundaryVector();
  resolveBottomBoundaryVector();
}

void FEMALEMethod::resolveLeftBoundaryVector() {
  const size_t celli = 0;
  switch (problem.leftBoundaryType) {
    case BoundaryType::eFree:
    case BoundaryType::eNoSlipWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = 0; i < (order + 1) * (order + 1); i += order + 1) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            const size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianX.coeffRef(indI, indJ) = 0.0;
            vectorLaplacianY.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = 0; i < (order + 1) * (order + 1); i += order + 1) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            const size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianX.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}

void FEMALEMethod::resolveTopBoundaryVector() {
  const size_t cellj = ySize - 1;
  switch (problem.topBoundaryType) {
    case BoundaryType::eFree:
    case BoundaryType::eNoSlipWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = order * (order + 1); i < (order + 1) * (order + 1);
             i++) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianX.coeffRef(indI, indJ) = 0.0;
            vectorLaplacianY.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = order * (order + 1); i < (order + 1) * (order + 1);
             i++) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianY.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}

void FEMALEMethod::resolveRightBoundaryVector() {
  const size_t celli = xSize - 1;
  switch (problem.rightBoundaryType) {
    case BoundaryType::eFree:
    case BoundaryType::eNoSlipWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = order; i < (order + 1) * (order + 1); i += order + 1) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianX.coeffRef(indI, indJ) = 0.0;
            vectorLaplacianY.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eWall:
      for (size_t cellj = 0; cellj < ySize; cellj++) {
        for (size_t i = order; i < (order + 1) * (order + 1); i += order + 1) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianX.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
  }
}

void FEMALEMethod::resolveBottomBoundaryVector() {
  const size_t cellj = 0;
  switch (problem.bottomBoundaryType) {
    case BoundaryType::eFree:
    case BoundaryType::eNoSlipWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = 0; i < order + 1; i++) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianX.coeffRef(indI, indJ) = 0.0;
            vectorLaplacianY.coeffRef(indI, indJ) = 0.0;
          }
        }
      }
      break;
    case BoundaryType::eWall:
      for (size_t celli = 0; celli < xSize; celli++) {
        for (size_t i = 0; i < order + 1; i++) {
          const size_t indI = getKinematicIndexFromCell(celli, cellj, i);
          for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
            size_t indJ = getKinematicIndexFromCell(celli, cellj, j);
            vectorLaplacianY.coeffRef(indI, indJ) = 0.0;
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

    Fe = Fx.transpose() * u05;
    Fe += Fy.transpose() * v05;
    e05 = Mt_inv * Fe;
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

  Fe = Fx.transpose() * u05;
  Fe += Fy.transpose() * v05;
  e05 = Mt_inv * Fe;
  e05 *= dt;
  e05 += e;

  x05 = x + dt * u05;
  y05 = y + dt * v05;

  t += dt;
  if (dt <= gamma * tau) {
    dt = beta2 * dt;
  }

  x.swap(x05);
  y.swap(y05);
  e.swap(e05);
}

void FEMALEMethod::remap() {
  optimizeMesh();
  transitionToRemap();
  remapT = remapTMin;
  while (remapT < remapTMax) {
    remapStep();
  }
  transitionToLagrange();
}

void FEMALEMethod::remapStep() {
  remapDt = remapTMax - remapT;

  calcRemapMatrices(x, y);
  if (remapTau < remapDt) {
    remapDt = remapTau;
  }
  u = MvSolver.solve(Kv * u);
  assert(MvSolver.info() == Eigen::Success);
  v = MvSolver.solve(Kv * v);
  assert(MvSolver.info() == Eigen::Success);

  remapStepRhoPrepare();
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      remapStepRhoMinMaxCell(celli, cellj);
      remapStepRhoAvgCell(celli, cellj);
    }
  }
  remapStepRhoFlux();
  remapStepRhoFinal();

  remapStepRhoEPrepare();
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      remapStepRhoEMinMaxCell(celli, cellj);
      remapStepRhoEAvgCell(celli, cellj);
    }
  }
  remapStepRhoEFlux();
  remapStepRhoEFinal();

  x += uMesh * remapDt;
  y += vMesh * remapDt;
}

void FEMALEMethod::remapStepRhoPrepare() {
  rhoRemapHigh = K * rhoRemap;
  rhoRemapHigh *= remapDt;
  rhoRemapHigh = MInv * rhoRemapHigh;
  rhoRemapHigh += rhoRemap;
  rhoRemapLow = Kupwinded * rhoRemap;
  rhoRemapLow *= remapDt;
  rhoRemapLow = MlumpedInv * rhoRemapLow;
  rhoRemapLow += rhoRemap;
}

void FEMALEMethod::remapStepRhoMinMaxCell(const size_t celli,
                                          const size_t cellj) {
  const double rhoMinCell = getMinRhoFromRemapCellAndNeighbours(celli, cellj);
  const double rhoMaxCell = getMaxRhoFromRemapCellAndNeighbours(celli, cellj);
  for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
    const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
    rhoRemapMin(indI) = rhoMinCell;
    rhoRemapMax(indI) = rhoMaxCell;
  }
}

void FEMALEMethod::remapStepRhoAvgCell(const size_t celli, const size_t cellj) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MCell(
      (order + 1) * (order + 1), (order + 1) * (order + 1));
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapLowCell((order + 1) *
                                                           (order + 1));

  for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
    const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
    rhoRemapLowCell(i) = rhoRemapLow(indI);
    for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
      const size_t indJ = getAdvectionIndexFromCell(celli, cellj, j);
      MCell(i, j) = M.coeff(indI, indJ);
    }
  }

  const double rhoAvgNumerator =
      Eigen::VectorXd::Ones((order + 1) * (order + 1)).transpose() * MCell *
      rhoRemapLowCell;
  const double rhoAvgDenominator =
      Eigen::VectorXd::Ones((order + 1) * (order + 1)).transpose() * MCell *
      Eigen::VectorXd::Ones((order + 1) * (order + 1));
  const double rhoAvgCell = rhoAvgNumerator / rhoAvgDenominator;

  for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
    const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
    rhoRemapAvg(indI) = rhoAvgCell;
  }
}

void FEMALEMethod::remapStepRhoFlux() {
  antidiffusiveFlux.setZero();
  antidiffusiveFlux.reserve(
      Eigen::VectorXi::Constant(Na, 4 * (order + 1) * (order + 1)));
  for (int k = 0; k < M.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      const double Mij = it.value();
      const double dRhoHighi = rhoRemapHigh(indI) - rhoRemap(indI);
      const double dRhoHighj = rhoRemapHigh(indJ) - rhoRemap(indJ);
      antidiffusiveFlux.coeffRef(indI, indJ) += Mij * (dRhoHighi - dRhoHighj);
    }
  }
  for (int k = 0; k < D.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      const double Dij = it.value();
      const double dRhoij = rhoRemap(indI) - rhoRemap(indJ);
      antidiffusiveFlux.coeffRef(indI, indJ) += remapDt * Dij * dRhoij;
    }
  }
  antidiffusiveFlux.makeCompressed();
  const size_t nInCell = (order + 1) * (order + 1);
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
        const double mi = Mlumped.coeff(i, i);
        const double dRhoLowi = rhoRemapLow(indI) - rhoRemapAvg(indI);
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          const size_t indJ = getAdvectionIndexFromCell(celli, cellj, j);
          const double mj = Mlumped.coeff(j, j);
          const double dRhoLowj = rhoRemapLow(indJ) - rhoRemapAvg(indJ);

          antidiffusiveFlux.coeffRef(indI, indJ) +=
              (mi * dRhoLowi - mj * dRhoLowj) / nInCell;
        }
      }
    }
  }
  assert(antidiffusiveFlux.isCompressed());
  remapStepRhoFluxLimiting();
}

void FEMALEMethod::remapStepRhoFluxLimiting() {
  fluxLimitingFactors.setZero();
  fluxLimitingFactors.reserve(
      Eigen::VectorXi::Constant(Na, 4 * (order + 1) * (order + 1)));
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoMax = rhoRemapMax - rhoRemapAvg;
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoMin = rhoRemapMin - rhoRemapAvg;
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoPlus(Na);
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoMinus(Na);
  for (int k = 0; k < antidiffusiveFlux.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(antidiffusiveFlux, k);
         it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      if (indI == indJ) {
        continue;
      }
      const double gij = it.value();
      if (gij >= 0.0) {
        dRhoPlus(indI) += gij / Mlumped.coeff(indI, indI);
      } else {
        dRhoMinus(indI) += gij / Mlumped.coeff(indI, indI);
      }
    }
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> betaPlus(Na);
  Eigen::Matrix<double, Eigen::Dynamic, 1> betaMinus(Na);
  for (size_t i = 0; i < Na; i++) {
    betaPlus(i) = std::min(1.0, dRhoMax(i) / dRhoPlus(i));
    betaMinus(i) = std::min(1.0, dRhoMin(i) / dRhoMinus(i));
  }

  for (int k = 0; k < antidiffusiveFlux.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(antidiffusiveFlux, k);
         it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      const double gij = it.value();
      if (gij > 0.0) {
        fluxLimitingFactors.coeffRef(indI, indJ) =
            std::min(betaPlus(indI), betaMinus(indJ));
      } else {
        fluxLimitingFactors.coeffRef(indI, indJ) =
            std::min(betaMinus(indI), betaPlus(indJ));
      }
    }
  }
  fluxLimitingFactors.makeCompressed();
}

void FEMALEMethod::remapStepRhoFinal() {
  rhoRemap.setZero();

  for (int k = 0; k < antidiffusiveFlux.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(antidiffusiveFlux, k);
         it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      if (indI == indJ) {
        continue;
      }
      const double gij = it.value();
      const double betaij = fluxLimitingFactors.coeff(indI, indJ);
      const double mi = Mlumped.coeff(indI, indI);
      rhoRemap(indI) += betaij * gij / mi;
    }
  }

  rhoRemap += rhoRemapAvg;
}

void FEMALEMethod::remapStepRhoEPrepare() {
  rhoERemapHigh = K * rhoERemap;
  rhoERemapHigh *= remapDt;
  rhoERemapHigh = MInv * rhoERemapHigh;
  rhoERemapHigh += rhoERemap;
  rhoERemapLow = Kupwinded * rhoERemap;
  rhoERemapLow *= remapDt;
  rhoERemapLow = MlumpedInv * rhoERemapLow;
  rhoERemapLow += rhoERemap;
}

void FEMALEMethod::remapStepRhoEMinMaxCell(const size_t celli,
                                           const size_t cellj) {
  const double eMinCell = getMinEFromRemapCellAndNeighbours(celli, cellj);
  const double eMaxCell = getMaxEFromRemapCellAndNeighbours(celli, cellj);
  for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
    const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
    eRemapMin(indI) = eMinCell;
    eRemapMin(indI) = eMaxCell;
  }
}

void FEMALEMethod::remapStepRhoEAvgCell(const size_t celli,
                                        const size_t cellj) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MCell(
      (order + 1) * (order + 1), (order + 1) * (order + 1));
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoRemapLowCell((order + 1) *
                                                           (order + 1));
  Eigen::Matrix<double, Eigen::Dynamic, 1> rhoERemapLowCell((order + 1) *
                                                            (order + 1));

  for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
    const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
    rhoRemapLowCell(i) = rhoRemapLow(indI);
    rhoERemapLowCell(i) = rhoERemapLow(indI);
    for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
      const size_t indJ = getAdvectionIndexFromCell(celli, cellj, j);
      MCell(i, j) = M.coeff(indI, indJ);
    }
  }

  const double eAvgNumerator =
      Eigen::VectorXd::Ones((order + 1) * (order + 1)).transpose() * MCell *
      rhoERemapLowCell;
  const double eAvgDenominator =
      Eigen::VectorXd::Ones((order + 1) * (order + 1)).transpose() * MCell *
      rhoRemapLowCell;
  const double eAvgCell = eAvgNumerator / eAvgDenominator;

  for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
    const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
    rhoERemapAvg(indI) = rhoRemap(indI) * eAvgCell;
  }
}

void FEMALEMethod::remapStepRhoEFlux() {
  antidiffusiveFlux.setZero();
  antidiffusiveFlux.reserve(
      Eigen::VectorXi::Constant(Na, 4 * (order + 1) * (order + 1)));
  for (int k = 0; k < M.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      const double Mij = it.value();
      const double dRhoEHighi = rhoERemapHigh(indI) - rhoERemap(indI);
      const double dRhoEHighj = rhoERemapHigh(indJ) - rhoERemap(indJ);
      antidiffusiveFlux.coeffRef(indI, indJ) += Mij * (dRhoEHighi - dRhoEHighj);
    }
  }
  for (int k = 0; k < D.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      const double Dij = it.value();
      const double dRhoEij = rhoERemap(indI) - rhoERemap(indJ);
      antidiffusiveFlux.coeffRef(indI, indJ) += remapDt * Dij * dRhoEij;
    }
  }
  antidiffusiveFlux.makeCompressed();
  const size_t nInCell = (order + 1) * (order + 1);
  for (size_t celli = 0; celli < xSize; celli++) {
    for (size_t cellj = 0; cellj < ySize; cellj++) {
      for (size_t i = 0; i < (order + 1) * (order + 1); i++) {
        const size_t indI = getAdvectionIndexFromCell(celli, cellj, i);
        const double mi = Mlumped.coeff(i, i);
        const double dRhoELowi = rhoERemapLow(indI) - rhoERemapAvg(indI);
        for (size_t j = 0; j < (order + 1) * (order + 1); j++) {
          const size_t indJ = getAdvectionIndexFromCell(celli, cellj, j);
          const double mj = Mlumped.coeff(j, j);
          const double dRhoELowj = rhoERemapLow(indJ) - rhoERemapAvg(indJ);

          antidiffusiveFlux.coeffRef(indI, indJ) +=
              (mi * dRhoELowi - mj * dRhoELowj) / nInCell;
        }
      }
    }
  }
  assert(antidiffusiveFlux.isCompressed());
  remapStepRhoEFluxLimiting();
}

void FEMALEMethod::remapStepRhoEFluxLimiting() {
  fluxLimitingFactors.setZero();
  fluxLimitingFactors.reserve(
      Eigen::VectorXi::Constant(Na, 4 * (order + 1) * (order + 1)));
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoEMax =
      (rhoRemap.array() * eRemapMax.array()).matrix() - rhoERemapAvg;
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoEMin =
      (rhoRemap.array() * eRemapMin.array()).matrix() - rhoERemapAvg;
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoEPlus(Na);
  Eigen::Matrix<double, Eigen::Dynamic, 1> dRhoEMinus(Na);
  for (int k = 0; k < antidiffusiveFlux.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(antidiffusiveFlux, k);
         it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      if (indI == indJ) {
        continue;
      }
      const double gij = it.value();
      if (gij >= 0.0) {
        dRhoEPlus(indI) += gij / Mlumped.coeff(indI, indI);
      } else {
        dRhoEMinus(indI) += gij / Mlumped.coeff(indI, indI);
      }
    }
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> betaPlus(Na);
  Eigen::Matrix<double, Eigen::Dynamic, 1> betaMinus(Na);
  for (size_t i = 0; i < Na; i++) {
    betaPlus(i) = std::min(1.0, dRhoEMax(i) / dRhoEPlus(i));
    betaMinus(i) = std::min(1.0, dRhoEMin(i) / dRhoEMinus(i));
  }

  for (int k = 0; k < antidiffusiveFlux.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(antidiffusiveFlux, k);
         it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      const double gij = it.value();
      if (gij > 0.0) {
        fluxLimitingFactors.coeffRef(indI, indJ) =
            std::min(betaPlus(indI), betaMinus(indJ));
      } else {
        fluxLimitingFactors.coeffRef(indI, indJ) =
            std::min(betaMinus(indI), betaPlus(indJ));
      }
    }
  }
  fluxLimitingFactors.makeCompressed();
}

void FEMALEMethod::remapStepRhoEFinal() {
  rhoERemap.setZero();

  for (int k = 0; k < antidiffusiveFlux.outerSize(); k++) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(antidiffusiveFlux, k);
         it; ++it) {
      const size_t indI = it.row();
      const size_t indJ = it.col();
      if (indI == indJ) {
        continue;
      }
      const double gij = it.value();
      const double betaij = fluxLimitingFactors.coeff(indI, indJ);
      const double mi = Mlumped.coeff(indI, indI);
      rhoERemap(indI) += betaij * gij / mi;
    }
  }

  rhoERemap += rhoERemapAvg;
}

Eigen::Matrix<double, 2, 2> FEMALEMethod::getCellJacobian(
    const size_t celli, const size_t cellj, const size_t quadi,
    const size_t quadj, const Eigen::Matrix<double, Eigen::Dynamic, 1> &x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y) const {
  Eigen::Matrix<double, 2, 2> jacobian = Eigen::Matrix<double, 2, 2>::Zero();
  for (size_t k = 0; k < (order + 1) * (order + 1); k++) {
    size_t indK = getKinematicIndexFromCell(celli, cellj, k);
    double xk = x(indK);
    double yk = y(indK);
    double basis2Ddx =
        kinematicBasis1DdxQuadValues[(k % (order + 1)) * quadOrder + quadi] *
        kinematicBasis1DQuadValues[(k / (order + 1)) * quadOrder + quadj];
    double basis2Ddy =
        kinematicBasis1DQuadValues[(k % (order + 1)) * quadOrder + quadi] *
        kinematicBasis1DdxQuadValues[(k / (order + 1)) * quadOrder + quadj];
    jacobian(0, 0) += xk * basis2Ddx;
    jacobian(0, 1) += xk * basis2Ddy;
    jacobian(1, 0) += yk * basis2Ddx;
    jacobian(1, 1) += yk * basis2Ddy;
  }
  return jacobian;
}

void FEMALEMethod::FCTProjection() {
  const size_t nProj = xProj.size();

  assert(wQuad.size() == quadOrder * quadOrder);
  assert(detJacQuad.size() == quadOrder * quadOrder);
  assert(basisQuadValues.rows() == nProj);
  assert(basisQuadValues.cols() == quadOrder * quadOrder);
  assert(xQuad.size() == quadOrder * quadOrder);
  assert(xProj.size() == nProj);
  assert(yQuad.size() == quadOrder * quadOrder);
  assert(ymin.size() == nProj);
  assert(ymax.size() == nProj);

  a.resize(nProj);
  a.setZero();
  b.resize(nProj);
  b.setZero();
  z.resize(nProj);
  z.setZero();
  FCTCoeffPlus.resize(nProj);
  FCTCoeffPlus.setZero();
  FCTCoeffMinus.resize(nProj);
  FCTCoeffMinus.setZero();
  localMassMatrix.resize(nProj, nProj);
  localMassMatrix.setZero();
  localMassMatrix_inv.resize(nProj, nProj);
  localMassMatrixLumped.resize(nProj, nProj);
  localMassMatrixLumped.reserve(Eigen::VectorXi::Constant(nProj, 1));
  localFluxMatrix.resize(nProj, nProj);
  FCTCoeffMatrix.resize(nProj, nProj);
  FCTCoeffMatrix.setZero();

  // Compute localMassMatrix, inverse, a, b
  for (size_t quad = 0; quad < quadOrder * quadOrder; quad++) {
    const double jacDet = detJacQuad(quad);
    const double weight = wQuad(quad);
    const double xi = xQuad(quad);
    const double yi = yQuad(quad);

    for (size_t basisi = 0; basisi < nProj; basisi++) {
      const double basisi2D = basisQuadValues(basisi, quad);
      for (size_t basisj = 0; basisj < nProj; basisj++) {
        const double basisj2D = basisQuadValues(basisj, quad);
        const double quad2D = weight * basisi2D * basisj2D * jacDet;
        localMassMatrix(basisi, basisj) += quad2D;
        localMassMatrixLumped.coeffRef(basisi, basisi) += quad2D;
      }
      a(basisi) += weight * xi * basisi2D * jacDet;
      b(basisi) += weight * xi * yi * basisi2D * jacDet;
    }
  }
  localMassMatrix_inv = localMassMatrix.inverse();
  localMassMatrixLumped.makeCompressed();

  const double yAvg = b.sum() / a.sum();
  xyLow = xProj * yAvg;
  xyHigh = localMassMatrix_inv * b;

  z = b;
  z -= localMassMatrixLumped * xyLow;

  for (size_t i = 0; i < nProj; i++) {
    for (size_t j = 0; j < nProj; j++) {
      localFluxMatrix(i, j) = localMassMatrix(i, j) * (xyHigh(i) - xyHigh(j)) +
                              (z(i) - z(j)) / nProj;
    }
  }

  // calc FCT coeffs
  for (size_t i = 0; i < nProj; i++) {
    const double deltaxyMax = xProj(i) * ymax(i) - xyLow(i);
    const double deltaxyMin = xProj(i) * ymin(i) - xyLow(i);
    double deltaxyPlus = 0.0;
    double deltaxyMinus = 0.0;
    for (size_t j = 0; j < nProj; j++) {
      if (j == i) {
        continue;
      }
      const double fij = localFluxMatrix(i, j);
      if (fij >= 0.0) {
        deltaxyPlus += fij;
      } else {
        deltaxyMinus += fij;
      }
    }
    deltaxyPlus /= localMassMatrixLumped.coeff(i, i);
    deltaxyMinus /= localMassMatrixLumped.coeff(i, i);

    FCTCoeffPlus(i) = std::min(1.0, deltaxyMax / deltaxyPlus);
    FCTCoeffMinus(i) = std::min(1.0, deltaxyMin / deltaxyMinus);
  }
  for (size_t i = 0; i < nProj; i++) {
    for (size_t j = 0; j < nProj; j++) {
      const double fij = localFluxMatrix(i, j);
      if (fij > 0.0) {
        FCTCoeffMatrix(i, j) = std::min(FCTCoeffPlus(i), FCTCoeffMinus(j));
      } else if (fij < 0.0) {
        FCTCoeffMatrix(i, j) = std::min(FCTCoeffMinus(i), FCTCoeffPlus(j));
      }
    }
  }

  xyOut = xyLow;
  for (size_t i = 0; i < nProj; i++) {
    double sum = 0.0;
    for (size_t j = 0; j < nProj; j++) {
      const double alphaij = FCTCoeffMatrix(i, j);
      const double fij = localFluxMatrix(i, j);

      sum += alphaij * fij;
    }
    xyOut(i) += sum;
  }
}
