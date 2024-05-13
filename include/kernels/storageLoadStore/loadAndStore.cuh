#pragma once
#include "basic_data/qcu_complex.cuh"
namespace qcu {
namespace kernel {
// COALESCED MEMORY
template <typename _Float>
static __device__ __forceinline__ void loadGauge(Complex<_Float> *u_local, void *gauge_ptr, int direction,
                                                 const Point<_Float> &p, int sub_Lx, int Ly, int Lz, int Lt, int Nc) {
  Complex<_Float> *u = p.getCoalescedGaugeAddr(gauge_ptr, direction, sub_Lx, Ly, Lz, Lt, Nc);
  int half_vol = sub_Lx * Ly * Lz * Lt;
  for (int i = 0; i < Nc * Nc; i++) {
    u_local[i] = *u;
    u += half_vol;
  }
}

// version 2 does not reconstruct the SU3 matrix, only stores the 12 complex numbers
template <typename _Float>
static __device__ __forceinline__ void loadGauge2(Complex<_Float> *u_local, void *gauge_ptr, int direction,
                                                  const Point<_Float> &p, int sub_Lx, int Ly, int Lz, int Lt, int Nc) {
  Complex<_Float> *u = p.getCoalescedGaugeAddr(gauge_ptr, direction, sub_Lx, Ly, Lz, Lt, Nc);
  int half_vol = sub_Lx * Ly * Lz * Lt;
  for (int i = 0; i < (Nc - 1) * Nc; i++) {
    u_local[i] = *u;
    u += half_vol;
  }
}

template <typename _Float>
static __device__ __forceinline__ void storeGauge(void *gauge_ptr, Complex<_Float> *u_local, int direction,
                                                  const Point<_Float> &p, int sub_Lx, int Ly, int Lz, int Lt, int Nc) {
  Complex<_Float> *u = p.getCoalescedGaugeAddr(gauge_ptr, direction, sub_Lx, Ly, Lz, Lt, Nc);
  int half_vol = sub_Lx * Ly * Lz * Lt;
  for (int i = 0; i < Nc * Nc; i++) {
    *u = u_local[i];
    u += half_vol;
  }
}
// loadVector<_Float>(srcLocal, fermionIn, move_point, Lx, Ly, Lz, Lt, _nColor);
// COALESCED
template <typename _Float>
static __device__ __forceinline__ void loadVector(Complex<_Float> *src_local, void *fermion_in, const Point<_Float> &p,
                                                  int sub_Lx, int Ly, int Lz, int Lt, int Nc) {
  Complex<_Float> *src = p.getCoalescedVectorAddr(fermion_in, sub_Lx, Ly, Lz, Lt);
  int half_vol = sub_Lx * Ly * Lz * Lt;
  for (int i = 0; i < Ns * Nc; i++) {
    src_local[i] = *src;
    src += half_vol;
  }
}

// COALESCED:
template <typename _Float>
static __device__ __forceinline__ void storeVector(Complex<_Float> *dst_local, void *fermion_out,
                                                   const Point<_Float> &p, int sub_Lx, int Ly, int Lz, int Lt, int Nc) {
  Complex<_Float> *src = p.getCoalescedVectorAddr(fermion_out, sub_Lx, Ly, Lz, Lt);
  int half_vol = sub_Lx * Ly * Lz * Lt;
  for (int i = 0; i < Ns * Nc; i++) {
    *src = dst_local[i];
    src += half_vol;
  }
}
}  // namespace kernel
}  // namespace qcu