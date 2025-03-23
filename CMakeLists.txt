cmake_minimum_required(VERSION 3.12)
project(MolecularSim)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files with wildcards
file(GLOB SRC_FILES src/mc/*.cpp)

# Include directories for header files
include_directories(src/mc)

# Add the executable
add_executable(mc_simulation ${SRC_FILES})

# Add the path to the data directory as a preprocessor definition
add_definitions(-DDATA_PATH="${CMAKE_SOURCE_DIR}/data")
