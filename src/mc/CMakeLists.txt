# Add source files for the `mc` version
file(GLOB MC_SRC_FILES *.cpp)

# Include headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Define executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(mc_simulation ${MC_SRC_FILES})

# Link libraries (if any) and preprocessor definitions
target_compile_definitions(mc_simulation PRIVATE DATA_PATH="${CMAKE_SOURCE_DIR}/data")
