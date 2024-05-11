#pragma once

#include "qcu_complex.cuh"
#include "qcu_enum.h"
namespace qcu {
namespace kernel {

// for su(N), donnot set nColor = 3
// _TComplex is means type of complex (_Complex is reserved for FC99)
template <template <typename> class _TComplex, typename _Float, int _dim, int _dir>
template <class _TComplex>
static __device__ __forceinline__ void dslashSpMV_4d_kernel(_TComplex<_Float> *u_local, _TComplex<_Float> *src_local,
                                                            _TComplex<_Float> *dst_local, _Float flag, int nColor,
                                                           ) {
  Complex temp1;
  Complex temp2;

  if (dim == X_DIM) {
    switch (dir) {
      case FWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1 = temp2 = 0;
#pragma unroll
          for (int j = 0; j < Nc; j++) {
            temp1 += (src_local[0 * Nc + j] - src_local[3 * Nc + j].multiply_i() * flag) * u_local[i * Nc + j];
            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] - src_local[2 * Nc + j].multiply_i() * flag) * u_local[i * Nc + j];
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[3 * Nc + i] += temp1.multiply_i() * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[2 * Nc + i] += temp2.multiply_i() * flag;
        }
      } break;
      case BWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();
#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] + src_local[3 * Nc + j].multiply_i() * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj

            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] + src_local[2 * Nc + j].multiply_i() * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[3 * Nc + i] += temp1.multiply_minus_i() * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[2 * Nc + i] += temp2.multiply_minus_i() * flag;
        }
      } break;
      default:
        break;
    }
  } else if (dim == Y_DIM) {
    switch (dir) {
      case FWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();

#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] + src_local[3 * Nc + j] * flag) * u_local[i * Nc + j];
            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] - src_local[2 * Nc + j] * flag) * u_local[i * Nc + j];
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[3 * Nc + i] += temp1 * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[2 * Nc + i] += -temp2 * flag;
        }
      } break;
      case BWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();
#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] - src_local[3 * Nc + j] * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
                                                  // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] + src_local[2 * Nc + j] * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[3 * Nc + i] += -temp1 * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[2 * Nc + i] += temp2 * flag;
        }
      } break;
      default:
        break;
    }
  } else if (dim == Z_DIM) {
    switch (dir) {
      case FWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();

#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] - src_local[2 * Nc + j].multiply_i() * flag) * u_local[i * Nc + j];
            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] + src_local[3 * Nc + j].multiply_i() * flag) * u_local[i * Nc + j];
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[2 * Nc + i] += temp1.multiply_i() * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[3 * Nc + i] += temp2.multiply_minus_i() * flag;
        }
      } break;
      case BWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();

#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] + src_local[2 * Nc + j].multiply_i() * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] - src_local[3 * Nc + j].multiply_i() * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[2 * Nc + i] += temp1.multiply_minus_i() * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[3 * Nc + i] += temp2.multiply_i() * flag;
        }
      } break;
      default:
        break;
    }
  } else if (dim == T_DIM) {
    switch (dir) {
      case FWD: {
#pragma unroll
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();

#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] - src_local[2 * Nc + j] * flag) * u_local[i * Nc + j];
            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] - src_local[3 * Nc + j] * flag) * u_local[i * Nc + j];
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[2 * Nc + i] += -temp1 * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[3 * Nc + i] += -temp2 * flag;
        }
      } break;
      case BWD: {
        for (int i = 0; i < Nc; i++) {
          temp1.clear2Zero();
          temp2.clear2Zero();

#pragma unroll
          for (int j = 0; j < Nc; j++) {
            // first row vector with col vector
            temp1 += (src_local[0 * Nc + j] + src_local[2 * Nc + j] * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
            // second row vector with col vector
            temp2 += (src_local[1 * Nc + j] + src_local[3 * Nc + j] * flag) *
                     u_local[j * Nc + i].conj();  // transpose and conj
          }
          dst_local[0 * Nc + i] += temp1;
          dst_local[2 * Nc + i] += temp1 * flag;
          dst_local[1 * Nc + i] += temp2;
          dst_local[3 * Nc + i] += temp2 * flag;
        }
      } break;
      default:
        break;
    }
  } else {
    assert(0);
  }
}

template <template <typename> class _TComplex, typename _Float>
__device__ __forceinline__ void dslashSpMV_4d_kernel<X_DIM, FWD>(_TComplex<_Float> *uLocal, _TComplex<_Float> *srcLocal,
                                                                 _TComplex<_Float> *dstLocal, _Float daggerFlag,
                                                                 int nColor) {
  Complex temp1;
  Complex temp2;

#pragma unroll
  for (int i = 0; i < nColor; i++) {
    temp1 = 0;
    temp2 = 0;

#pragma unroll
    for (int j = 0; j < nColor; j++) {
      temp1 += (srcLocal[0 * nColor + j] - srcLocal[3 * nColor + j].multiply_i() * daggerFlag) * uLocal[i * nColor + j];
      // second row vector with col vector
      temp2 += (srcLocal[1 * nColor + j] - srcLocal[2 * nColor + j].multiply_i() * daggerFlag) * uLocal[i * nColor + j];
    }
    dstLocal[0 * nColor + i] += temp1;
    dstLocal[3 * nColor + i] += temp1.multiply_i() * daggerFlag;
    dstLocal[1 * nColor + i] += temp2;
    dstLocal[2 * nColor + i] += temp2.multiply_i() * daggerFlag;
  }
}

template <template <typename> class _TComplex, typename _Float>
__device__ __forceinline__ void dslashSpMV_4d_kernel<X_DIM, BWD>(_TComplex<_Float> *uLocal, _TComplex<_Float> *srcLocal,
                                                                 _TComplex<_Float> *dstLocal, _Float daggerFlag,
                                                                 int nColor) {
  Complex temp1;
  Complex temp2;

#pragma unroll
  for (int i = 0; i < nColor; i++) {
    temp1 = 0;
    temp2 = 0;
#pragma unroll
    for (int j = 0; j < nColor; j++) {
      // first row vector with col vector
      temp1 += (srcLocal[0 * nColor + j] + srcLocal[3 * nColor + j].multiply_i() * daggerFlag) *
               uLocal[j * nColor + i].conj();  // transpose and conj

      // second row vector with col vector
      temp2 += (srcLocal[1 * nColor + j] + srcLocal[2 * nColor + j].multiply_i() * daggerFlag) *
               uLocal[j * nColor + i].conj();  // transpose and conj
    }
    dstLocal[0 * nColor + i] += temp1;
    dstLocal[3 * nColor + i] += temp1.multiply_minus_i() * daggerFlag;
    dstLocal[1 * nColor + i] += temp2;
    dstLocal[2 * nColor + i] += temp2.multiply_minus_i() * daggerFlag;
  }
}

}  // namespace kernel

}  // namespace qcu