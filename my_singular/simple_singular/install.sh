#!/bin/bash

# Exit script on error
set -e

# Create necessary directories
mkdir -p ~/mySingular/my_singular/simple_singular/build_dir ~/mySingular//my_singular/simple_singular/install_dir

# Define variables
INSTALL_PREFIX="/home/santosh/mySingular/my_singular/simple_singular/install_dir"
BUILD_TYPE="Release"
BOOST_NO_CMAKE="on"
BUILD_DIR="/home/santosh/mySingular/my_singular/simple_singular/build_dir"
SOURCE_DIR="/home/santosh/mySingular/my_singular/simple_singular/"

# Set Singular installation directory
SINGULAR_INSTALL_DIR="/home/santosh/singular-gpispace/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/singular-4.4.0p2-syrkttc4im2j3tzob5jykruuxnushksj"

# Export the environment variable for CMake to use
export SINGULAR_INSTALL_DIR
export LD_LIBRARY_PATH="$SINGULAR_INSTALL_DIR/lib:$LD_LIBRARY_PATH"



# Remove old build directory and create a fresh one
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

# Run CMake
cmake -D CMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
      -D CMAKE_BUILD_TYPE="$BUILD_TYPE" \
      -B "$BUILD_DIR" \
      -S "$SOURCE_DIR"

# Build and install
cmake --build "$BUILD_DIR" --target install -- -j "$(nproc)"
