cmake_minimum_required(VERSION 3.14)
project(gravity_core)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG     v2.11.1
)
FetchContent_MakeAvailable(pybind11)

#Note: we name the module "gravity_core" to match the C++ code
pybind11_add_module(gravity_core src/engine.cpp)

if(MSVC)
    target_compile_options(gravity_core PRIVATE /O2)
else()
    target_compile_options(gravity_core PRIVATE -O3 -march=native)
endif()
