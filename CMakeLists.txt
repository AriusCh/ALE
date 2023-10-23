cmake_minimum_required(VERSION 3.17..3.27)

project(ale-solver VERSION 0.1
    DESCRIPTION "An Arbitrary Lagrangian-Euleran solver"
    LANGUAGES C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(SOURCE_FILES src/main.cpp
    src/method.hpp src/method.cpp
    src/primitive.hpp src/primitive.cpp
)

add_executable(ale ${SOURCE_FILES})

set_target_properties(ale PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQIRED YES
    CXX_EXTENSIONS NO)


