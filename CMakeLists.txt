cmake_minimum_required(VERSION 3.12)
project(MolecularSim)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include subdirectories
add_subdirectory(src/mc)
add_subdirectory(src/mcMPI)
add_subdirectory(src/mcCU)
