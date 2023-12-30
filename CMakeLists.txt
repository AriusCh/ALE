cmake_minimum_required(VERSION 3.17..3.28)

project(ale-solver VERSION 0.1
    DESCRIPTION "An Arbitrary Lagrangian-Euleran solver"
    LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(SOURCE_FILES src/main.cpp
    src/simulation.hpp src/simulation.cpp
    src/method.hpp src/method.cpp
    src/grid.hpp src/grid.cpp
    src/boundary.hpp src/boundary.cpp
    src/problem.hpp src/problem.cpp
    src/output_mgr.hpp src/output_mgr.cpp
    src/eos.hpp src/eos.cpp
    src/logger.hpp src/logger.cpp
)

add_executable(ale ${SOURCE_FILES})

set_target_properties(ale PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQIRED YES
    CXX_EXTENSIONS NO)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
target_link_libraries(ale PUBLIC MPI::MPI_CXX)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
target_link_libraries(ale PUBLIC Eigen3::Eigen)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(ale PUBLIC OpenMP::OpenMP_CXX)
endif()

# if (MSVC) 
#     target_compile_options(ale PRIVATE /W4 /WX)
# else()
target_compile_options(ale PRIVATE -Wall -Wextra -Wpedantic -Werror)
# endif()

