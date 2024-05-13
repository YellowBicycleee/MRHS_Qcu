#pragma once

#include <assert.h>

#include "basic_data/qcu_complex.cuh"
#include "qcu.h"
#include "qcu_enum.h"
#include "qcu_macro.cuh"
#define QCU_PERFORMANCE
// _Complex keyword reserved for C99, donnot use it

namespace qcu {

template <typename _Float>
class Point {
 private:
  int x_;
  int y_;
  int z_;
  int t_;
  int parity_;

 public:
  Point() = default;
  __device__ __forceinline__ Point(const Point &rhs)
      : x_(rhs.x_), y_(rhs.y_), z_(rhs.z_), t_(rhs.t_), parity_(rhs.parity_) {}

  __device__ __forceinline__ Point(int x, int y, int z, int t, int parity)
      : x_(x), y_(y), z_(z), t_(t), parity_(parity) {}

  __device__ __forceinline__ int getParity() const { return parity_; }

  __device__ __forceinline__ Point move(int front_back, int direction, int Lx, int Ly, int Lz,
                                        int Lt) const {  // direction +-1234
#ifndef QCU_PERFORMANCE
    assert(direction == X_DIM || direction == Y_DIM || direction == Z_DIM || direction == T_DIM);
    assert(front_back == BWD || front_back == FWD);
#endif

    int new_pos;
    int eo = (y_ + z_ + t_) & 0x01;  // (y+z+t)%2

    if (direction == X_DIM) {
      if (!front_back) {  // front_back == BWD
        new_pos = x_ + (eo == parity_) * (-1 + (x_ == 0) * Lx);
        return Point(new_pos, y_, z_, t_, 1 - parity_);
      } else {  // front_back == FWD
        new_pos = x_ + (eo != parity_) * (1 + (x_ == Lx - 1) * (-Lx));
        return Point(new_pos, y_, z_, t_, 1 - parity_);
      }
    } else if (direction == Y_DIM) {
      if (!front_back) {  // front_back == BWD
        new_pos = y_ - 1 + (y_ == 0) * Ly;
        return Point(x_, new_pos, z_, t_, 1 - parity_);
      } else {  // front_back == FWD
        new_pos = y_ + 1 + (y_ == Ly - 1) * (-Ly);
        return Point(x_, new_pos, z_, t_, 1 - parity_);
      }
    } else if (direction == Z_DIM) {
      if (!front_back) {
        new_pos = z_ - 1 + (z_ == 0) * Lz;
        return Point(x_, y_, new_pos, t_, 1 - parity_);
      } else {
        new_pos = z_ + 1 + (z_ == Lz - 1) * (-Lz);
        return Point(x_, y_, new_pos, t_, 1 - parity_);
      }
    } else if (direction == T_DIM) {
      if (!front_back) {
        new_pos = t_ - 1 + (t_ == 0) * Lt;
        return Point(x_, y_, z_, new_pos, 1 - parity_);
      } else {
        new_pos = t_ + 1 + (t_ == Lt - 1) * (-Lt);
        return Point(x_, y_, z_, new_pos, 1 - parity_);
      }
    } else {
      return *this;
    }
  }

  __device__ __forceinline__ Complex<_Float> *getCoalescedVectorAddr(void *origin, int Lx, int Ly, int Lz,
                                                                     int Lt) const {
    return static_cast<Complex<_Float> *>(origin) + (((t_ * Lz + z_) * Ly + y_) * Lx + x_);
  }
  __device__ __forceinline__ Complex<_Float> *getCoalescedGaugeAddr(void *origin, int direction, int sub_Lx, int Ly,
                                                                    int Lz, int Lt, int Nc) const {
    return static_cast<Complex<_Float> *>(origin) + (direction * 2 + parity_) * sub_Lx * Ly * Lz * Lt * Nc * Nc +
           (((t_ * Lz + z_) * Ly + y_) * sub_Lx + x_);
  }

  __device__ __forceinline__ Point &operator=(const Point &rhs) {
    x_ = rhs.x_;
    y_ = rhs.y_;
    z_ = rhs.z_;
    t_ = rhs.t_;
    parity_ = rhs.parity_;
    return *this;
  }
  __device__ __forceinline__ Complex<_Float> *getPointClover(Complex<_Float> *origin, int Lx, int Ly, int Lz, int Lt,
                                                             int Nc) const {
    return origin + (((((parity_ * Lt + t_) * Lz + z_) * Ly + y_) * Lx) + x_) * (Nc * Ns * Nc * Ns / 2);
  }

  /// @brief get the vector from ptr to dst
  /// @param dst the destination vector, contiuously stored
  /// @param ptr the source vector, not contiuously stored, the stride is the distance between two elements
  /// @param elementNums the number of elements in the vector
  /// @param startPos the start position of the vector
  /// @param stride the distance between two elements
  __device__ __forceinline__ void getVector(Complex<_Float> *dst, Complex<_Float> *ptr, int elementNums, int startPos,
                                            int stride) {
    for (int i = 0; i < elementNums; i++) {
      dst[i] = ptr[startPos + i * stride];
    }
  }
};

}  // namespace qcu