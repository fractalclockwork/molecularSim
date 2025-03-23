# Enable CUDA support
enable_language(CUDA)

# Add source files for mcCU
file(GLOB MCCU_SRC_FILES *.cpp *.cu)

# Include directories for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Define CUDA executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(mdCU_simulation ${MCCU_SRC_FILES})

# Add nvcc compile options
set_target_properties(mdCU_simulation PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_ARCHITECTURES 60 # Replace "60" with your target architecture
    CUDA_NVCC_FLAGS "--fmad=false" # Think about this.
)

# Link necessary libraries (if any)
target_link_libraries(mdCU_simulation m)

# Link libraries (if any) and preprocessor definitions
target_compile_definitions(mdCU_simulation PRIVATE DATA_PATH="${CMAKE_SOURCE_DIR}/data")
