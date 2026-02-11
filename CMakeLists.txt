cmake_minimum_required(VERSION 3.14)
project(gravity_core)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(OpenMP REQUIRED)

# Physics Library
add_library(physics_lib INTERFACE)
target_include_directories(physics_lib INTERFACE src)
target_compile_options(physics_lib INTERFACE -O3 -march=native)
target_link_libraries(physics_lib INTERFACE OpenMP::OpenMP_CXX)

# Python Bindings
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(gravity_core src/engine.cpp)
target_link_libraries(gravity_core PRIVATE physics_lib)

# Testing
enable_testing()
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.9
)
FetchContent_MakeAvailable(Catch2)

if(EXISTS "${catch2_SOURCE_DIR}/contrib")
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")
endif()

add_executable(unit_tests tests/test_physics.cpp)
target_link_libraries(unit_tests PRIVATE physics_lib Catch2::Catch2)

include(CTest)
include(Catch)
catch_discover_tests(unit_tests)
