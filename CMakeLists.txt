cmake_minimum_required(VERSION 3.8)
include_directories(/usr/lib/x86_64-linux-gnu/openmpi/include)
set(CMAKE_CUDA_ARCHITECTURES 70)    # for V100
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_C_STANDARD 11)
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -fsanitize=address")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g")
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-gencode;arch=compute_70,code=sm_70")



project(QCU LANGUAGES C CXX CUDA) # project name = ${PROJECT_NAME}

add_library(qcu 
    SHARED
    ${PROJECT_SOURCE_DIR}/src/storage/qcu_storage.cu
    ${PROJECT_SOURCE_DIR}/src/qcd/qcu_wilson.cu   
    ${PROJECT_SOURCE_DIR}/src/qcu_dirac.cu
    ${PROJECT_SOURCE_DIR}/src/qcu.cu
)
add_library(qcu1
    STATIC 
    ${PROJECT_SOURCE_DIR}/src/storage/qcu_storage.cu
    ${PROJECT_SOURCE_DIR}/src/qcd/qcu_wilson.cu   
    ${PROJECT_SOURCE_DIR}/src/qcu_dirac.cu
    ${PROJECT_SOURCE_DIR}/src/qcu.cu
)

target_link_libraries(qcu
	cudart
	cudadevrt
)
target_link_libraries(qcu1
	cudart
	cudadevrt
)
target_include_directories(qcu 
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_include_directories(qcu1 
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
add_executable(test
    ${PROJECT_SOURCE_DIR}/src/test.cu
    # ${PROJECT_SOURCE_DIR}/tests/test_arg/test_arg.cu
)
target_include_directories(test 
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(test 
	qcu1
#	-lcudart
#	-lcudadevrt
)
# set_target_properties(test_linear_algebra PROPERTIES
    # RUNTIME_OUTPUT_DIRECTORY "tests")

# find_package(MPI REQUIRED)
# target_link_libraries(qcu ${MPI_LIBRARIES})
# target_link_libraries(qcu nccl)
