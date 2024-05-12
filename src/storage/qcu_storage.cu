#include "kernels/storage/storage_kernels.cuh"
#include "qcu_enum.h"
#include "storage/qcu_storage.h"

// TODO: 模板特化（现在只有SU3， double）
namespace qcu {
namespace storage {

void shiftVectorStorage(void *dstVec, void *srcVec, ShiftDirection shift_direction, int subLx, int Ly, int Lz, int Lt) {
  int half_vol = subLx * Ly * Lz * Lt;

  int block_size = BLOCK_SIZE;
  int grid_size = (half_vol + block_size - 1) / block_size;

  if (shift_direction == TO_COALESCE) {
    kernel::storeVectorCoalesced<double, 3><<<grid_size, block_size>>>(dstVec, srcVec, subLx, Ly, Lz, Lt);
  } else {
    kernel::storeVectorNonCoalesced<double, 3><<<grid_size, block_size>>>(dstVec, srcVec, subLx, Ly, Lz, Lt);
  }
  CHECK_CUDA(cudaDeviceSynchronize());
}

void shiftGaugeStorage(void *dst_vec, void *src_vec, ShiftDirection shift_direction, int Lx, int Ly, int Lz, int Lt) {
  int vol = Lx * Ly * Lz * Lt;
  int half_vol = vol / 2;
  int block_size = BLOCK_SIZE;
  int grid_size = (half_vol + block_size - 1) / block_size;
  if (shift_direction == TO_COALESCE) {
    kernel::storeGaugeCoalesced<double, 3><<<grid_size, block_size>>>(dst_vec, src_vec, Lx, Ly, Lz, Lt);
    CHECK_CUDA(cudaDeviceSynchronize());
  }
}
}  // namespace storage
}  // namespace qcu