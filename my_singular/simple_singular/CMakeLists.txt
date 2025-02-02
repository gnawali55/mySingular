cmake_minimum_required(VERSION 3.10)
project(simple_singular VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Manually set Singular paths
set(SINGULAR_INSTALL_DIR "$ENV{SINGULAR_INSTALL_DIR}")

set(SINGULAR_INCLUDE_DIRS 
    "${SINGULAR_INSTALL_DIR}/include"
    "${SINGULAR_INSTALL_DIR}/include/singular"
    "${SINGULAR_INSTALL_DIR}/include/singular/Singular"
)

set(SINGULAR_LIBRARY_DIR "${SINGULAR_INSTALL_DIR}/lib")

# Debugging - print paths
message(STATUS "SINGULAR_INSTALL_DIR: ${SINGULAR_INSTALL_DIR}")
message(STATUS "SINGULAR_INCLUDE_DIRS: ${SINGULAR_INCLUDE_DIRS}")
message(STATUS "SINGULAR_LIBRARY_DIR: ${SINGULAR_LIBRARY_DIR}")

# Include Singular headers
include_directories(${SINGULAR_INCLUDE_DIRS})

# Add executable
add_executable(simple_singular src/main.cpp)

# Link all Singular libraries
file(GLOB SINGULAR_LIBRARIES "${SINGULAR_LIBRARY_DIR}/*.so")
target_link_libraries(simple_singular ${SINGULAR_LIBRARIES})
target_link_libraries(simple_singular -lpthread -ldl -lm)

# Install the executable
install(TARGETS simple_singular
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)