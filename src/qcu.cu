#include <cuda.h>

#include "qcu.h"
#include "qcuDesc/qcu_desc.h"
#include "qcu_dirac.h"
#include "qcu_enum.h"
#include "qcu_kernel_param.cuh"
namespace qcu {

static int nColor = 3;  // number of colors, SU(N)
static ExecutionStreams gpuStreams;

class Qcu {
  //  protected:
  //   DSLASH_TYPE dslashType_;
  QcuLattDesc<4> lattDesc_;
  QcuProcDesc<4> procDesc_;
  //   Dslash *dslash;  // dslash operator
  Dirac *dirac;

 public:
  Qcu(QcuParam *pLattDesc, QcuGrid *pProcDesc) : lattDesc_(pLattDesc), procDesc_(pProcDesc), dirac(nullptr) {}
  void getDslash(DSLASH_TYPE dslashType, void *gauge, double mass, int nColors, int nInputs,
                 QCU_PRECISION floatPrecision) {
    switch (dslashType) {
      case DSLASH_WILSON: {
        dirac = new DiracWilson(gauge, mass, nColors, nInputs, floatPrecision, lattDesc_, procDesc_, dslashType,
                                gpuStreams);
      } break;
      default:
        errorQcu("Now we only support Wilson dslash.\n");
        break;
    }
  }
  void qcuDslash(void* outputMRHS, void* inputMRHS, int parity, int daggerFlag) {
    dirac->Dslash(outputMRHS, inputMRHS, parity, daggerFlag);
  }

  ~Qcu() {
    if (dirac) {
      delete dirac;
      dirac = nullptr;
    }
  }
};

}  // namespace qcu

static qcu::Qcu *qcu_ptr = nullptr;

void initGridSize(QcuGrid_t *grid, QcuParam *p_param, void *gauge, void *fermion_in, void *fermion_out) {}

void destroyQcu() {}

void dslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int parity) {}

void fullDslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int dagger_flag) {}
void cg_inverter(void *x_vector, void *b_vector, void *gauge, QcuParam *param, double p_max_prec, double p_kappa) {}

void loadQcuGauge(void *gauge, QcuParam *param) {}