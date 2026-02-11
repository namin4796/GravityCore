cmake_minimum_required(VERSION 3.14)

# Force macOS settings to be modern
if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum OS X deployment version" FORCE)
endif()

project(gravity_core)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. OpenMP Setup
find_package(OpenMP REQUIRED)

# 2. Physics Library
add_library(physics_lib INTERFACE)
target_include_directories(physics_lib INTERFACE src)
target_compile_options(physics_lib INTERFACE -O3 -march=native)
target_link_libraries(physics_lib INTERFACE OpenMP::OpenMP_CXX)

# 3. Python Module
include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(gravity_core src/engine.cpp)
target_link_libraries(gravity_core PRIVATE physics_lib)

# 4. Catch2 (The New Testing Framework)
enable_testing()

# We download the single-header version directly. No compilation needed.
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.9
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/contrib")

# 5. Test Executable
add_executable(unit_tests tests/test_physics.cpp)
target_link_libraries(unit_tests PRIVATE physics_lib Catch2::Catch2)

# Enable CTest integration
include(CTest)
include(Catch)
catch_discover_tests(unit_tests)
