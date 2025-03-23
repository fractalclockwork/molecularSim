#Add source files for the `mcMPI` version
file(GLOB MCMPI_SRC_FILES*.cpp)

#Include headers
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Define executable
        set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR} / bin)
            add_executable(mcMPI_simulation ${MCMPI_SRC_FILES})

#Link MPI library
                find_package(MPI REQUIRED) include_directories(${MPI_INCLUDE_PATH})
                    target_link_libraries(mcMPI_simulation ${MPI_LIBRARIES})

#Link additional libraries
                        target_link_libraries(mcMPI_simulation m)

#Preprocessor definition for data path
                            target_compile_definitions(
                                mcMPI_simulation PRIVATE DATA_PATH = "${CMAKE_SOURCE_DIR}/data")
