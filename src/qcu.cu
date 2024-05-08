#include <cuda.h>

#include "qcu.h"
#include "qcuDesc/qcuDesc.h"

namespace qcu {

static int nColor = 3;  // number of colors, SU(N)

class Qcu {
 protected:
  QcuLattDesc<4> lattDesc;
  QcuProcDesc<4> procDesc;

 public:
  Qcu(QcuGrid_t *pLattDesc, QcuParam *pProcDesc) : lattDesc(pLattDesc), procDesc(pProcDesc) {}
  ~Qcu() {}
};

}  // namespace qcu

static qcu::Qcu *qcu_ptr = nullptr;

void initGridSize(QcuGrid_t *grid, QcuParam *p_param, void *gauge, void *fermion_in, void *fermion_out) {}

void destroyQcu() {}

void dslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int parity) {}

void fullDslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int dagger_flag) {}
void cg_inverter(void *x_vector, void *b_vector, void *gauge, QcuParam *param, double p_max_prec, double p_kappa) {}

void loadQcuGauge(void *gauge, QcuParam *param) {}