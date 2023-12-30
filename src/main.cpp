#include <mpi.h>

#include <Eigen/Core>
#include <format>
#include <iostream>
#include <memory>
#include <string>

#include "logger.hpp"
#include "method.hpp"
#include "problem.hpp"
#include "simulation.hpp"

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  MPI_Init(&argc, &argv);

  int worldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  char processorName[MPI_MAX_PROCESSOR_NAME];
  int processorNameLen;
  MPI_Get_processor_name(processorName, &processorNameLen);

  Logger mainLogger;
  std::string initMessage =
      std::format("NODE: {}, PROCESS: {}, WORLD_SIZE: {}, OPENMP THREADS: {}",
                  processorName, worldRank, worldSize, Eigen::nbThreads());
  mainLogger.log(initMessage, LogLevel::eGeneral);

  MPI_Finalize();
  return 0;
}
