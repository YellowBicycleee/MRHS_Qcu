#pragma once
#include "qcu_enum.h"
namespace qcu {
namespace storage {
void shiftVectorStorage(void *dstVec, void *srcVec, ShiftDirection shift_direction, int subLx, int Ly, int Lz, int Lt);
void shiftGaugeStorage(void *dst_vec, void *src_vec, ShiftDirection shift_direction, int Lx, int Ly, int Lz, int Lt);
}  // namespace storage
}  // namespace qcu