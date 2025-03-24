#Add source files for the `mc` version
file(GLOB MC_SRC_FILES *.cpp)

#Include headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Define executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(mc_simulation ${MC_SRC_FILES})

# Fetch Toml++
include(FetchContent)
FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(tomlplusplus)

message("Toml++ source directory: ${tomlplusplus_SOURCE_DIR}")

# Add the tomlplusplus headers to your target
target_include_directories(mc_simulation PRIVATE ${tomlplusplus_SOURCE_DIR}/include)

#Link libraries(if any) and preprocessor definitions
target_compile_definitions(mc_simulation PRIVATE DATA_PATH="${CMAKE_SOURCE_DIR}/data")
