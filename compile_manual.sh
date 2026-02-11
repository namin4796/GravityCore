#!/bin/bash

# Clean
rm -rf build
mkdir -p build

# Get Python configuration
PY_INCLUDES=$(python3 -m pybind11 --includes)
PY_EXT_SUFFIX=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
OMP_PREFIX=$(brew --prefix libomp)

echo "Compiling GravityCore manually..."

# Compile
clang++ -O3 -Wall -shared -std=c++17 -undefined dynamic_lookup \
    $PY_INCLUDES \
    -I src \
    -Xpreprocessor -fopenmp \
    -I"$OMP_PREFIX/include" \
    -L"$OMP_PREFIX/lib" -lomp \
    src/engine.cpp \
    -o "build/gravity_core$PY_EXT_SUFFIX"

if [ $? -eq 0 ]; then
    echo "Build Successful!"
else
    echo "Build Failed."
    exit 1
fi
