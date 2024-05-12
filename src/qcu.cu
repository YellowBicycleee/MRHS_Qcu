#include <cuda.h>

#include <cstdio>
#include <vector>

#include "qcu.h"
#include "qcuDesc/qcu_desc.h"
#include "qcu_dirac.h"
#include "qcu_enum.h"
#include "qcu_kernel_param.cuh"
#include "storage/qcu_storage.h"
using std::vector;
namespace qcu {

static ExecutionStreams gpuStreams;
class Qcu {
  //  protected:
  //   DSLASH_TYPE dslashType_;
  QcuLattDesc<4> lattDesc_;
  QcuProcDesc<4> procDesc_;
  //   Dslash *dslash;  // dslash operator
  Dirac *dirac_;
  vector<void *> fermionIn_queue;
  vector<void *> fermionOut_queue;
  int nColors_;
  void *localGauge_;

 public:
  Qcu(QcuParam *pLattDesc, QcuGrid *pProcDesc, int nColors)
      : lattDesc_(pLattDesc), procDesc_(pProcDesc), dirac_(nullptr), localGauge_(nullptr), nColors_(nColors) {}
  void getDslash(DSLASH_TYPE dslashType, void *gauge, double mass, int nColors, int nInputs,
                 QCU_PRECISION floatPrecision) {
    switch (dslashType) {
      case DSLASH_WILSON: {
        dirac_ = new DiracWilson(gauge, mass, nColors, floatPrecision, lattDesc_, procDesc_, dslashType, gpuStreams);
        // DiracWilson(gauge, mass, nColors, floatPrecision, lattDesc, procDesc, dslashType, streams)
      } break;
      default:
        errorQcu("Now we only support Wilson dslash.\n");
        break;
    }
  }

  void startDslash(int parity, int daggerFlag = 0) {
    // copy the input and output fermion pointers to the corresponding arrays
    // void *outputMRHS_ptr = static_cast<void *>(fermionOut_queue.data());
    // void *inputMRHS_ptr = static_cast<void *>(fermionIn_queue.data());

    vector<void *> coalescedOutputMRHS(fermionOut_queue.size());
    vector<void *> coalescedInputMRHS(fermionIn_queue.size());
    for (int i = 0; i < fermionOut_queue.size(); i++) {
      CHECK_CUDA(cudaMalloc(&coalescedOutputMRHS[i], 2 * sizeof(double) * lattDesc_.X() / 2 * lattDesc_.Y() *
                                                         lattDesc_.Z() * lattDesc_.T() * Ns * nColors_));
      CHECK_CUDA(cudaMalloc(&coalescedInputMRHS[i], 2 * sizeof(double) * lattDesc_.X() / 2 * lattDesc_.Y() *
                                                        lattDesc_.Z() * lattDesc_.T() * Ns * nColors_));
      storage::shiftVectorStorage(coalescedInputMRHS[i], fermionIn_queue[i], TO_COALESCE, lattDesc_.X() / 2,
                                  lattDesc_.Y(), lattDesc_.Z(), lattDesc_.T());
    }
    void *outputMRHS_ptr = static_cast<void *>(coalescedOutputMRHS.data());
    void *inputMRHS_ptr = static_cast<void *>(coalescedInputMRHS.data());

    dirac_->Dslash(outputMRHS_ptr, inputMRHS_ptr, parity, daggerFlag, fermionOut_queue.size());

    for (int i = 0; i < fermionIn_queue.size(); i++) {
      storage::shiftVectorStorage(fermionOut_queue[i], coalescedOutputMRHS[i], TO_NON_COALESCE, lattDesc_.X() / 2,
                                  lattDesc_.Y(), lattDesc_.Z(), lattDesc_.T());
      CHECK_CUDA(cudaFree(coalescedOutputMRHS[i]));
      CHECK_CUDA(cudaFree(coalescedInputMRHS[i]));
    }
    fermionIn_queue.clear();
    fermionOut_queue.clear();
  }
  void pushFermion(void *fermionOut, void *fermionIn) {
    fermionIn_queue.push_back(fermionIn);
    fermionOut_queue.push_back(fermionOut);
  }
  void loadGauge(void *gauge) {
    if (localGauge_ == nullptr) {
      //   errorQcu("Gauge is already loaded.\n");
      CHECK_CUDA(cudaMalloc(&localGauge_, 2 * sizeof(double) * lattDesc_.X() * lattDesc_.Y() * lattDesc_.Z() *
                                              lattDesc_.T() * Nd * nColors_ * nColors_));
    }
    storage::shiftGaugeStorage(localGauge_, gauge, TO_COALESCE, lattDesc_.X(), lattDesc_.Y(), lattDesc_.Z(),
                               lattDesc_.T());
  }

  ~Qcu() {
    if (dirac_) {
      CHECK_CUDA(cudaFree(localGauge_));
      delete dirac_;
      dirac_ = nullptr;
    }
  }
};

}  // namespace qcu

static qcu::Qcu *qcu_ptr = nullptr;

void initGridSize(QcuGrid_t *grid, QcuParam *p_param, int nColors) {
  if (qcu_ptr != nullptr) {
    delete qcu_ptr;
    qcu_ptr = nullptr;

    qcu_ptr = new qcu::Qcu(p_param, grid, nColors);
  }
}

void pushBackFermions(void *fermionOut, void *fermionIn) {  // 押入输入输出向量
  if (qcu_ptr == nullptr) {
    errorQcu("Qcu is not initialized.\n");
  }
  qcu_ptr->pushFermion(fermionOut, fermionIn);
}

void startDslash(int parity, int daggerFlag) {
  if (qcu_ptr == nullptr) {
    errorQcu("Qcu is not initialized.\n");
  }
  qcu_ptr->startDslash(parity, daggerFlag);
}

void loadQcuGauge(void *gauge) {
  if (qcu_ptr == nullptr) {
    errorQcu("Qcu is not initialized.\n");
  }
  qcu_ptr->loadGauge(gauge);
}
// void loadQcuGauge(void *gauge, QcuParam *param); // 废弃

void dslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int parity) {}           // 废弃
void fullDslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int dagger_flag) {}  // 废弃
// TODO
void cg_inverter(void *x_vector, void *b_vector, void *gauge, QcuParam *param, double p_max_prec, double p_kappa) {}

__attribute__((destructor)) void destroyQcu() {
  if (qcu_ptr != nullptr) {
    delete qcu_ptr;
    qcu_ptr = nullptr;
  }
}