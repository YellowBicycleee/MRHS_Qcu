#pragma once

namespace qcu {

template <size_t _Nvec>
struct MrhsInput {
  void* fermion[_Nvec];
  void* data() { return fermion; }
  __device__ __host__ void* operator[](size_t i) { return fermion[i]; }
};

};  // namespace qcu