#pragma once

// #include <cuda.h>
#include <cuda_runtime.h>
#include "qcu_macro.cuh"
#include <cstdio>
namespace qcu {

struct ExecutionStreams {
  cudaStream_t stream;
  cudaStream_t boundaryStreams[4 * 2];  // 4 dims, 2 dirs

  ExecutionStreams() {
    CHECK_CUDA(cudaStreamCreate(&stream));
    for (int i = 0; i < 4 * 2; i++) {
      CHECK_CUDA(cudaStreamCreate(&boundaryStreams[i]));
    }
  }
  ExecutionStreams(const ExecutionStreams &execStreams) {
    stream = execStreams.stream;
    for (int i = 0; i < 4 * 2; i++) {
      boundaryStreams[i] = execStreams.boundaryStreams[i];
    }
  }
  void operator=(const ExecutionStreams &execStreams) {
    stream = execStreams.stream;
    for (int i = 0; i < 4 * 2; i++) {
      boundaryStreams[i] = execStreams.boundaryStreams[i];
    }
  }

  ~ExecutionStreams() {
    CHECK_CUDA(cudaStreamDestroy(stream));
    for (int i = 0; i < 4 * 2; i++) {
      CHECK_CUDA(cudaStreamDestroy(boundaryStreams[i]));
    }
  }
};

// struct BaseParam {
//   static ExecutionStreams gpuStreams;
// };

// ExecutionStreams BaseParam::gpuStreams;

// struct KernelParam : public BaseParam {
//   //   int blockSize[3];
//   //   int gridSize[3];
//   KernelParam() {}
//   //   {
//   // for (int i = 0; i < 3; i++) {
//   //   blockSize[i] = 1;
//   //   gridSize[i] = 1;
//   // }
//   // blockSize[0] = 256;
//   //   }
//   void setParam(int blockSize, int gridSize) {}
//   //   {  // 只用blockDim.x和gridDim.x
//   // blockSize[0] = x;
//   // gridSize[0] = y;
// //   //   }
// //   dim3 getBlockSize() const { return dim3(blockSize[0]); }  //{ return dim3(blockSize[0], blockSize[1], blockSize[2]); }
// //   dim3 getGridSize() const {}   //{ return dim3(gridSize[0], gridSize[1], gridSize[2]); }
// };

}  // namespace qcu