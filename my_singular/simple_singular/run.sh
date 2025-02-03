#!/bin/bash

# Set the Singular installation directory
SINGULAR_INSTALL_DIR="/home/santosh/singular-gpispace/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/singular-4.4.0p2-syrkttc4im2j3tzob5jykruuxnushksj"

# Set the library path
export LD_LIBRARY_PATH="$SINGULAR_INSTALL_DIR/lib:$LD_LIBRARY_PATH"

# Export the Singular installation directory for use in other commands
export SINGULAR_INSTALL_DIR

# Navigate to the directory containing the built executable
cd /home/santosh/mySingular/my_singular/simple_singular/install_dir/bin

# Run the executable with Valgrind for memory checking
echo "Running the program with Valgrind..."
valgrind --main-stacksize=16384 --tool=memcheck --leak-check=full --track-origins=yes --verbose ./simple_singular
# valgrind  --leak-check=full --show-leak-kinds=all ./simple_singular
gdb ./simple_singular
run
