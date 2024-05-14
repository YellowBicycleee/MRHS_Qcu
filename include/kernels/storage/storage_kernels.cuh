#pragma once

#include <cstdio>

#include "basic_data/qcu_complex.cuh"
#include "qcu_enum.h"
#include "qcu_macro.cuh"
#define WARP_SIZE 32
#define BLOCK_SIZE 256

namespace qcu {
namespace kernel {

template <typename _Float, int _nColor>
__global__ void storeVectorCoalesced(void* dstVec, void* srcVec, int subLx, int Ly, int Lz, int Lt) {
  constexpr int Nc = _nColor;

  int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int sub_vol = subLx * Ly * Lz * Lt;

  Complex<_Float>* dstVecPtr;
  Complex<_Float>* srcVecPtr;

  srcVecPtr = static_cast<Complex<_Float>*>(srcVec) + thread_id * Ns * Nc;
  dstVecPtr = static_cast<Complex<_Float>*>(dstVec) + thread_id;

  for (int i = 0; i < Ns * Nc; i++) {
    *dstVecPtr = *srcVecPtr;
    dstVecPtr += sub_vol;
    srcVecPtr++;
  }
}

template <typename _Float, int _nColor>
__global__ void storeVectorNonCoalesced(void* dstVec, void* srcVec, int subLx, int Ly, int Lz, int Lt) {
  constexpr int Nc = _nColor;

  int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int sub_vol = subLx * Ly * Lz * Lt;

  Complex<_Float>* dstVecPtr;
  Complex<_Float>* srcVecPtr;

  dstVecPtr = static_cast<Complex<_Float>*>(dstVec) + thread_id * Ns * Nc;
  srcVecPtr = static_cast<Complex<_Float>*>(srcVec) + thread_id;

  for (int i = 0; i < Ns * Nc; i++) {
    *dstVecPtr = *srcVecPtr;
    dstVecPtr++;
    srcVecPtr += sub_vol;
  }
}

template <typename _Float, int _nColor>
__global__ void storeGaugeCoalesced(void* dst_gauge, void* src_gauge, int Lx, int Ly, int Lz, int Lt) {
  constexpr int Nc = _nColor;

  int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int sub_Lx = Lx >> 1;
  int sub_vol = sub_Lx * Ly * Lz * Lt;

  Complex<_Float>* dstGaugePtr;
  Complex<_Float>* srcGaugePtr;

  for (int i = 0; i < Nd; i++) {
    for (int parity = 0; parity < 2; parity++) {
      dstGaugePtr = static_cast<Complex<_Float>*>(dst_gauge) + (2 * i + parity) * sub_vol * Nc * Nc + thread_id;
      srcGaugePtr = static_cast<Complex<_Float>*>(src_gauge) + (2 * i + parity) * sub_vol * Nc * Nc + thread_id * Nc * Nc;
      for (int j = 0; j < Nc * Nc; j++) {
        *dstGaugePtr = *srcGaugePtr;
        dstGaugePtr += sub_vol;
        srcGaugePtr++;
      }
    }
  }
}

template <typename _Float, int _nColor>
__global__ void storeGaugeNonCoalesced(void* dst_gauge, void* src_gauge, int Lx, int Ly, int Lz, int Lt) {
  constexpr int Nc = _nColor;

  int thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int sub_Lx = Lx >> 1;
  int sub_vol = sub_Lx * Ly * Lz * Lt;

  Complex<_Float>* dstGaugePtr;
  Complex<_Float>* srcGaugePtr;

  for (int i = 0; i < Nd; i++) {
    for (int parity = 0; parity < 2; parity++) {
      srcGaugePtr = static_cast<Complex<_Float>*>(src_gauge) + (2 * i + parity) * sub_vol * Nc * Nc + thread_id;
      dstGaugePtr =
          static_cast<Complex<_Float>*>(dst_gauge) + (2 * i + parity) * sub_vol * Nc * Nc + thread_id * Nc * Nc;
      for (int j = 0; j < Nc * Nc; j++) {
        *dstGaugePtr = *srcGaugePtr;
        srcGaugePtr += sub_vol;
        dstGaugePtr++;
      }
    }
  }
}

}  // namespace kernel
}  // namespace qcu