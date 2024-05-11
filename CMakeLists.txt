cmake_minimum_required(VERSION 3.8)
include_directories(/usr/lib/x86_64-linux-gnu/openmpi/include)
set(CMAKE_CUDA_ARCHITECTURES 70)    # for V100
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_C_STANDARD 11)

project(QCU LANGUAGES C CXX CUDA) # project name = ${PROJECT_NAME}

add_library(qcu 
    SHARED
    ${PROJECT_SOURCE_DIR}/src/qcu.cu
    ${PROJECT_SOURCE_DIR}/src/qcd/qcu_wilson.cu   
    ${PROJECT_SOURCE_DIR}/src/qcu_dirac.cu
)

target_include_directories(qcu 
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

# add_executable(test_linear_algebra 
    # ${PROJECT_SOURCE_DIR}/src/tests/linear_algebra/test_kernels.cu
    # ${PROJECT_SOURCE_DIR}/tests/test_arg/test_arg.cu
# )
# target_link_libraries(test_linear_algebra qcu)
# set_target_properties(test_linear_algebra PROPERTIES
    # RUNTIME_OUTPUT_DIRECTORY "tests")

# find_package(MPI REQUIRED)
# target_link_libraries(qcu ${MPI_LIBRARIES})
# target_link_libraries(qcu nccl)