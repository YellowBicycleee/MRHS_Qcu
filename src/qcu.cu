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
#include "timer.h"

// #define DEBUG

namespace qcu {

static ExecutionStreams gpuStreams;
class Qcu {
  //  protected:
  //   DSLASH_TYPE dslashType_;
  QcuLattDesc<4> lattDesc_;
  QcuProcDesc<4> procDesc_;
  //   Dslash *dslash;  // dslash operator
  Dirac *dirac_;
  vector<void *> fermionIn_queue_;
  vector<void *> fermionOut_queue_;

  vector<void *> coalesced_fermionIn_queue_;
  vector<void *> coalesced_fermionOut_queue_;

  int nColors_;
  void *localGauge_;

 public:
  Qcu(QcuParam *pLattDesc, QcuGrid *pProcDesc, int nColors)
      : lattDesc_(pLattDesc), procDesc_(pProcDesc), dirac_(nullptr), localGauge_(nullptr), nColors_(nColors) {}
  // TODO: 废弃gauge
  // void getDslash(DSLASH_TYPE dslashType, void *gauge, double mass, int nColors, int nInputs,
  //  QCU_PRECISION floatPrecision) {
  void getDslash(DSLASH_TYPE dslashType, double mass, int nColors, int nInputs, QCU_PRECISION floatPrecision) {
    if (localGauge_ == nullptr) {
      errorQcu("localGauge_ is nullptr, not loaded.\n");
    }
    switch (dslashType) {
      case DSLASH_WILSON: {
        dirac_ =
            new DiracWilson(localGauge_, mass, nColors, floatPrecision, lattDesc_, procDesc_, dslashType, gpuStreams);
        // DiracWilson(gauge, mass, nColors, floatPrecision, lattDesc, procDesc, dslashType, streams)
      } break;
      default:
        errorQcu("Now we only support Wilson dslash.\n");
        break;
    }
  }

  void startDslash(int parity, int daggerFlag = 0) {
    // copy the input and output fermion pointers to the corresponding arrays
    int inputNum = coalesced_fermionIn_queue_.size();
    void *outputMRHS_ptr = static_cast<void *>(coalesced_fermionOut_queue_.data());
    void *inputMRHS_ptr = static_cast<void *>(coalesced_fermionIn_queue_.data());

    dirac_->Dslash(outputMRHS_ptr, inputMRHS_ptr, parity, daggerFlag, fermionOut_queue_.size());

    for (int i = 0; i < inputNum; i++) {
      storage::shiftVectorStorage(fermionOut_queue_[i], coalesced_fermionOut_queue_[i], TO_NON_COALESCE,
                                  lattDesc_.X() / 2, lattDesc_.Y(), lattDesc_.Z(), lattDesc_.T());
    }
    for (int i = 0; i < coalesced_fermionIn_queue_.size(); i++) {
      void *fIn = coalesced_fermionIn_queue_[i];
      void *fOut = coalesced_fermionOut_queue_[i];
      CHECK_CUDA(cudaFree(fIn));
      CHECK_CUDA(cudaFree(fOut));
      coalesced_fermionIn_queue_[i] = nullptr;
      coalesced_fermionOut_queue_[i] = nullptr;
    }

#ifdef DEBUG
    printf("memory shift Over\n");
#endif
    fermionIn_queue_.clear();
    fermionOut_queue_.clear();
    coalesced_fermionIn_queue_.clear();
    coalesced_fermionOut_queue_.clear();
  }

  void pushFermion(void *fermionOut, void *fermionIn) {
    int halfVol = lattDesc_.X() / 2 * lattDesc_.Y() * lattDesc_.Z() * lattDesc_.T();
    int halfVectorLength = halfVol * Ns * nColors_;
    void *fIn;
    void *fOut;
    CHECK_CUDA(cudaMalloc(&fIn, 2 * sizeof(double) * halfVectorLength));
    CHECK_CUDA(cudaMalloc(&fOut, 2 * sizeof(double) * halfVectorLength));
    storage::shiftVectorStorage(fIn, fermionIn, TO_COALESCE, lattDesc_.X() / 2, lattDesc_.Y(), lattDesc_.Z(),
                                lattDesc_.T());
    coalesced_fermionIn_queue_.push_back(fIn);
    coalesced_fermionOut_queue_.push_back(fOut);
    fermionIn_queue_.push_back(fermionIn);
    fermionOut_queue_.push_back(fermionOut);
  }
  void loadGauge(void *gauge) {
    if (localGauge_ == nullptr) {
      CHECK_CUDA(cudaMalloc(&localGauge_, 2 * sizeof(double) * lattDesc_.X() * lattDesc_.Y() * lattDesc_.Z() *
                                              lattDesc_.T() * Nd * nColors_ * nColors_));
    }
    storage::shiftGaugeStorage(localGauge_, gauge, TO_COALESCE, lattDesc_.X(), lattDesc_.Y(), lattDesc_.Z(),
                               lattDesc_.T());
  }

  ~Qcu() {
    if (dirac_) {
      delete dirac_;
      dirac_ = nullptr;
    }
    if (localGauge_ != nullptr) {
      CHECK_CUDA(cudaFree(localGauge_));
      localGauge_ = nullptr;
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
  } else {
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

// void getDslash(int dslashType, void *gauge, double mass, int nColors, int nInputs, int floatPrecision) {
void getDslash(int dslashType, double mass, int nColors, int nInputs, int floatPrecision) {
  if (qcu_ptr == nullptr) {
    errorQcu("Qcu is not initialized.\n");
  }
  // qcu_ptr->getDslash(static_cast<DSLASH_TYPE>(dslashType), gauge, mass, nColors, nInputs,
  //                    static_cast<QCU_PRECISION>(floatPrecision));
  qcu_ptr->getDslash(static_cast<DSLASH_TYPE>(dslashType), mass, nColors, nInputs,
                     static_cast<QCU_PRECISION>(floatPrecision));
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
