#pragma once

#include <cassert>

#include "basic_data/qcu_complex.cuh"
#include "basic_data/qcu_point.cuh"
#include "kernels/storageLoadStore/loadAndStore.cuh"
#include "qcuDesc/mrhs_desc.h"
// #include "wrapper/type_wrapper.h"

namespace qcu {
namespace kernel {

// inline __device__ void calcCoord(int& x, int& y, int& z, int& t, int id, int Lx, int Ly, int Lz, int Lt) {
//   t = tid / (Lz * Ly * Lx);
//   z = tid % (Lz * Ly * Lx) / (Ly * Lx);
//   y = tid % (Ly * Lx) / Lx;
//   x = tid % Lx;
// }

// 命名规则：
// 模板参数以下划线开头，using 别名也用下划线开头
// class成员变量 以下划线结尾
// constexpr int Ndim = 4;
template <typename _Float, int _nColor, int _nInput>
// flag is daggerFlag
static __device__ void applyWilson(MrhsInput<_nInput>& outPtrs, MrhsInput<_nInput>& inPtrs, void* gauge,
                                   Complex<_Float>* regDst, Complex<_Float>* regSrc, Complex<_Float>* regGauge,
                                   _Float flag, int Lx, int Ly, int Lz, int Lt, int grid_x, int grid_y, int grid_z,
                                   int grid_t, int parity, int tid) {
  using Complex = Complex<_Float>;
  constexpr int Nc = _nColor;

  int coordBoundary;  // check if it is boundary of the lattice space
  int subLx = Lx >> 1;
  int t = tid / (Lz * Ly * subLx);
  int z = tid % (Lz * Ly * subLx) / (Ly * subLx);
  int y = tid % (Ly * subLx) / subLx;
  int x = tid % subLx;

  Point<_Float> p(x, y, z, t, parity);
  Point<_Float> move_point;
  int eo = (y + z + t) & 0x01;

  Complex temp1;
  Complex temp2;

  void* fermionIn;
  void* fermionOut;

  Complex* srcLocal = static_cast<Complex*>(regSrc);
  Complex* dstLocal;
  Complex* uLocal = static_cast<Complex*>(regGauge);
#pragma unroll
  for (int j = 0; j < Ns * _nColor * _nInput; j++) {
    regDst[j] = 0;
  }

  // X_DIM, FWD
  coordBoundary = (grid_x > 1 && x == subLx - 1 && parity != eo) ? subLx - 1 : subLx;
  loadGauge<_Float>(uLocal, gauge, X_DIM, p, subLx, Ly, Lz, Lt, _nColor);
  move_point = p.move(FWD, X_DIM, subLx, Ly, Lz, Lt);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (x < coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = temp2 = 0;

#pragma unroll
        for (int j = 0; j < Nc; j++) {
          temp1 += (srcLocal[0 * Nc + j] - srcLocal[3 * Nc + j].multiply_i() * flag) * uLocal[i * Nc + j];
          // second row vector with col vector
          temp2 += (srcLocal[1 * Nc + j] - srcLocal[2 * Nc + j].multiply_i() * flag) * uLocal[i * Nc + j];
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[3 * Nc + i] += temp1.multiply_i() * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[2 * Nc + i] += temp2.multiply_i() * flag;
      }
    }
  }

//   // X_DIM, BWD
  coordBoundary = (grid_x > 1 && x == 0 && parity == eo) ? 1 : 0;
  move_point = p.move(BWD, X_DIM, subLx, Ly, Lz, Lt);
  loadGauge(uLocal, gauge, X_DIM, move_point, subLx, Ly, Lz, Lt, _nColor);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (x >= coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;
#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 += (srcLocal[0 * Nc + j] + srcLocal[3 * Nc + j].multiply_i() * flag) *
                   uLocal[j * Nc + i].conj();  // transpose and conj

          // second row vector with col vector
          temp2 += (srcLocal[1 * Nc + j] + srcLocal[2 * Nc + j].multiply_i() * flag) *
                   uLocal[j * Nc + i].conj();  // transpose and conj
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[3 * Nc + i] += temp1.multiply_minus_i() * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[2 * Nc + i] += temp2.multiply_minus_i() * flag;
      }
    }
  }

  // Y_DIM  FWD
  coordBoundary = (grid_y > 1) ? Ly - 1 : Ly;
  loadGauge(uLocal, gauge, Y_DIM, p, subLx, Ly, Lz, Lt, _nColor);
  move_point = p.move(FWD, Y_DIM, subLx, Ly, Lz, Lt);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (y < coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;

#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 += (srcLocal[0 * Nc + j] + srcLocal[3 * Nc + j] * flag) * uLocal[i * Nc + j];
          // second row vector with col vector
          temp2 += (srcLocal[1 * Nc + j] - srcLocal[2 * Nc + j] * flag) * uLocal[i * Nc + j];
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[3 * Nc + i] += temp1 * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[2 * Nc + i] += -temp2 * flag;
      }
    }
  }

  // Y_DIM  BWD
  coordBoundary = (grid_y > 1) ? 1 : 0;
  move_point = p.move(BWD, Y_DIM, subLx, Ly, Lz, Lt);
  loadGauge(uLocal, gauge, Y_DIM, move_point, subLx, Ly, Lz, Lt, _nColor);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (y >= coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;
#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 += (srcLocal[0 * Nc + j] - srcLocal[3 * Nc + j] * flag) *
                   uLocal[j * Nc + i].conj();  // transpose and conj
                                               // second row vector with col vector
          temp2 +=
              (srcLocal[1 * Nc + j] + srcLocal[2 * Nc + j] * flag) * uLocal[j * Nc + i].conj();  // transpose and conj
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[3 * Nc + i] += -temp1 * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[2 * Nc + i] += temp2 * flag;
      }
    }
  }

  // Z_DIM  FWD
  coordBoundary = (grid_z > 1) ? Lz - 1 : Lz;
  loadGauge(uLocal, gauge, Z_DIM, p, subLx, Ly, Lz, Lt, _nColor);
  move_point = p.move(FWD, Z_DIM, subLx, Ly, Lz, Lt);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (z < coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;

#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 += (srcLocal[0 * Nc + j] - srcLocal[2 * Nc + j].multiply_i() * flag) * uLocal[i * Nc + j];
          // second row vector with col vector
          temp2 += (srcLocal[1 * Nc + j] + srcLocal[3 * Nc + j].multiply_i() * flag) * uLocal[i * Nc + j];
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[2 * Nc + i] += temp1.multiply_i() * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[3 * Nc + i] += temp2.multiply_minus_i() * flag;
      }
    }
  }

  // Z_DIM BWD
  coordBoundary = (grid_z > 1) ? 1 : 0;
  move_point = p.move(BWD, Z_DIM, subLx, Ly, Lz, Lt);
  loadGauge(uLocal, gauge, Z_DIM, move_point, subLx, Ly, Lz, Lt, _nColor);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (z >= coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;

#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 += (srcLocal[0 * Nc + j] + srcLocal[2 * Nc + j].multiply_i() * flag) *
                   uLocal[j * Nc + i].conj();  // transpose and conj
          // second row vector with col vector
          temp2 += (srcLocal[1 * Nc + j] - srcLocal[3 * Nc + j].multiply_i() * flag) *
                   uLocal[j * Nc + i].conj();  // transpose and conj
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[2 * Nc + i] += temp1.multiply_minus_i() * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[3 * Nc + i] += temp2.multiply_i() * flag;
      }
    }
  }

  // T_DIM  FWD
  coordBoundary = (grid_t > 1) ? Lt - 1 : Lt;
  loadGauge(uLocal, gauge, T_DIM, p, subLx, Ly, Lz, Lt, _nColor);
  move_point = p.move(FWD, T_DIM, subLx, Ly, Lz, Lt);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (t < coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;

#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 += (srcLocal[0 * Nc + j] - srcLocal[2 * Nc + j] * flag) * uLocal[i * Nc + j];
          // second row vector with col vector
          temp2 += (srcLocal[1 * Nc + j] - srcLocal[3 * Nc + j] * flag) * uLocal[i * Nc + j];
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[2 * Nc + i] += -temp1 * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[3 * Nc + i] += -temp2 * flag;
      }
    }
  }

  // T_DIM  BWD
  coordBoundary = (grid_t > 1) ? 1 : 0;
  move_point = p.move(BWD, T_DIM, subLx, Ly, Lz, Lt);
  loadGauge(uLocal, gauge, T_DIM, move_point, subLx, Ly, Lz, Lt, _nColor);
  for (int _ = 0; _ < _nInput; _++) {
    fermionIn = inPtrs[_];
    dstLocal = regDst + _ * Ns * Nc;
    loadVector<_Float>(srcLocal, fermionIn, move_point, subLx, Ly, Lz, Lt, _nColor);
    if (t >= coordBoundary) {
#pragma unroll
      for (int i = 0; i < Nc; i++) {
        temp1 = 0;
        temp2 = 0;
#pragma unroll
        for (int j = 0; j < Nc; j++) {
          // first row vector with col vector
          temp1 +=
              (srcLocal[0 * Nc + j] + srcLocal[2 * Nc + j] * flag) * uLocal[j * Nc + i].conj();  // transpose and conj
          // second row vector with col vector
          temp2 +=
              (srcLocal[1 * Nc + j] + srcLocal[3 * Nc + j] * flag) * uLocal[j * Nc + i].conj();  // transpose and conj
        }
        dstLocal[0 * Nc + i] += temp1;
        dstLocal[2 * Nc + i] += temp1 * flag;
        dstLocal[1 * Nc + i] += temp2;
        dstLocal[3 * Nc + i] += temp2 * flag;
      }
    }
  }

  // store
  for (int i = 0; i < _nInput; i++) {
    fermionOut = outPtrs[i];
    dstLocal = regDst + i * Ns * Nc;
    storeVector<_Float>(dstLocal, fermionOut, p, subLx, Ly, Lz, Lt, _nColor);
  }
}

// template <typename _Float>
template <typename _Float, int _nColor, int _nInput>
static void __global__ dslashKernel(MrhsInput<_nInput> outPtrs, MrhsInput<_nInput> inPtrs, void* gauge,
                                    _Float daggerFlag, int Lx, int Ly, int Lz, int Lt, int grid_x, int grid_y,
                                    int grid_z, int grid_t, int parity) {
  using Complex = qcu::Complex<_Float>;
  constexpr int Nc = _nColor;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = gridDim.x * blockDim.x;
  int lattVol = Lx * Ly * Lz * Lt / 2;
  Complex localSrc[Ns * Nc];
  Complex localDst[Ns * Nc * _nInput];
  Complex localGauge[Nc * Nc];

  for (int i = 0; i < lattVol; i += stride) {
    applyWilson<_Float, _nColor, _nInput>(outPtrs, inPtrs, gauge, localDst, localSrc, localGauge, daggerFlag, Lx, Ly,
                                          Lz, Lt, grid_x, grid_y, grid_z, grid_t, parity, tid);
  }
}

}  // namespace kernel

}  // namespace qcu