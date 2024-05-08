#pragma once
#include <basic_data/qcu_complex.cuh>

#include "qcu_enum.h"
namespace qcu {
namespace wrapper {

template <QCU_DAGGER_FLAG _Flag, typename _Float>
struct DaggerFlag {
  static constexpr _Float daggerFlag = (_Float)0;
};

template <typename _Float>
struct DaggerFlag<QCU_DAGGER_NO, _Float> {
  static constexpr _Float daggerFlag = (_Float)1.0;
};

template <typename _Float>
struct DaggerFlag<QCU_DAGGER_YES, _Float> {
  static constexpr _Float daggerFlag = (_Float)-1.0;
};

// float precision
template <QCU_PRECESION _Precision>
struct ComplexTypeWrapper {
  using Float = float;
  using Complex = qcu::Complex<float>;
};
template <>
struct ComplexTypeWrapper<QCU_HALF_PRECISION> {
  using Float = half;
  using Complex = qcu::Complex<half>;
};
template <>
struct ComplexTypeWrapper<QCU_DOUBLE_PRECISION> {
  using Float = double;
  using Complex = qcu::Complex<double>;
};

}  // namespace wrapper
}  // namespace qcu