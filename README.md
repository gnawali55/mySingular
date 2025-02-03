This is purely for Singular C++ code for computing Schyerer's frame. 
<!-- install.sh withb valgrind -->
#!/bin/bash

# Exit script on error
set -e

# Create necessary directories
mkdir -p ~/mySingular/my_singular/simple_singular/build_dir ~/mySingular/my_singular/simple_singular/install_dir

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

# Add the optional step to run with valgrind
read -p "Do you want to run the program with valgrind for memory checking? (y/n): " RUN_VALGRIND

if [ "$RUN_VALGRIND" == "y" ]; then
    echo "Running program with valgrind..."
    
    # Set path to your program here. If your binary is `simple_singular`, adjust accordingly.
    # This assumes the binary is in the install directory after building.
    PROGRAM_PATH="$INSTALL_PREFIX/bin/simple_singular"
    
    if [ -f "$PROGRAM_PATH" ]; then
        # Run your program with valgrind
        valgrind --tool=memcheck --leak-check=full --track-origins=yes --verbose "$PROGRAM_PATH"
    else
        echo "Error: $PROGRAM_PATH does not exist. Make sure the program was built successfully."
        exit 1
    fi
else
    echo "Skipping valgrind run."
fi
