#pragma once


// constexpr int Nc = 3;    // Nc isn't defined to 3
constexpr int Ns = 4;

#define CHECK_MPI(cmd)                         \
  do {                                         \
    int err = cmd;                             \
    if (err != MPI_SUCCESS) {                  \
      fprintf(stderr, "MPI error: %d\n", err); \
      exit(1);                                 \
    }                                          \
  } while (0)

#define CHECK_CUDA(cmd)                                                                                   \
  do {                                                                                                    \
    cudaError_t err = cmd;                                                                                \
    if (err != cudaSuccess) {                                                                             \
      fprintf(stderr, "CUDA error: %s, file %s, line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
      exit(1);                                                                                            \
    }                                                                                                     \
  } while (0)

#define CHECK_NCCL(cmd)                                             \
  do {                                                              \
    ncclResult_t err = cmd;                                         \
    if (err != ncclSuccess) {                                       \
      fprintf(stderr, "NCCL error: %s\n", ncclGetErrorString(err)); \
      exit(1);                                                      \
    }                                                               \
  } while (0)


#define errorQcu(msg) \
  do {                \
    fprintf(stderr, msg); \
    fprintf(stderr, "Error happened in file %s, line %d\n", __FILE__, __LINE__); \
    exit(1);          \
  } while (0)
