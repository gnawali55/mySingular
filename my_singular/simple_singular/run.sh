#!/bin/bash

# Set the Singular installation directory
SINGULAR_INSTALL_DIR="/home/santosh/singular-gpispace/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/singular-4.4.0p2-syrkttc4im2j3tzob5jykruuxnushksj"

# Set the library path
export LD_LIBRARY_PATH="$SINGULAR_INSTALL_DIR/lib:$LD_LIBRARY_PATH"

# Export the Singular installation directory for use in other commands
export SINGULAR_INSTALL_DIR

# Navigate to the directory containing the built executable
cd /home/santosh/mySingular/my_singular/simple_singular/install_dir/bin

# Run the executable

./simple_singular
# gdb ./simple_singular
# run
# bt  # Get backtrace when it crashes

