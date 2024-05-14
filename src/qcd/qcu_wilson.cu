#include <cuda_fp16.h>

#include "kernels/wilson_dslash_4d.cuh"
#include "qcd/qcu_dslash.h"
namespace qcu {

// using namespace qcu::kernel;

template <typename _Float, int _nColor, int _nInput>
static inline void executeWilson(MrhsInput<_nInput>& out, MrhsInput<_nInput>& in, void* gauge, _Float daggerFlag,
                                 int Lx, int Ly, int Lz, int Lt, int grid_x, int grid_y, int grid_z, int grid_t,
                                 int parity, int blockSize, cudaStream_t stream) {
  int gridSize = (Lx * Ly * Lz * Lt / 2 + blockSize - 1) / blockSize;
  kernel::dslashKernel<_Float, _nColor, _nInput><<<gridSize, blockSize, 0, stream>>>(out, in, gauge, daggerFlag, Lx, Ly, Lz, Lt,
                                                                             grid_x, grid_y, grid_z, grid_t, parity);
  CHECK_CUDA(cudaGetLastError()); 
  CHECK_CUDA(cudaDeviceSynchronize());
}

// TODO: 改blockSize不写死
template <typename _Float, int _nColor, int _nInput>
void WilsonDslash<_Float, _nColor, _nInput>::apply() {
  // int nInput = _nInput;
  int nColor = _nColor;
  int blockSize = 256;

  MrhsInput<_nInput>& out = *static_cast<MrhsInput<_nInput>*>(this->dslashParam_.fermionOut_MRHS);
  MrhsInput<_nInput>& in = *static_cast<MrhsInput<_nInput>*>(this->dslashParam_.fermionIn_MRHS);
  int Lx = this->dslashParam_.lattDesc.dim[0];
  int Ly = this->dslashParam_.lattDesc.dim[1];
  int Lz = this->dslashParam_.lattDesc.dim[2];
  int Lt = this->dslashParam_.lattDesc.dim[3];
  int grid_x = this->dslashParam_.procDesc.dim[0];
  int grid_y = this->dslashParam_.procDesc.dim[1];
  int grid_z = this->dslashParam_.procDesc.dim[2];
  int grid_t = this->dslashParam_.procDesc.dim[3];

  int parity = this->dslashParam_.parity;
  void* gauge = this->dslashParam_.gauge;

  _Float daggerFlag = (_Float)((this->dslashParam_.daggerFlag == 0) ? 1.0 : -1.0);

  executeWilson<_Float, _nColor, _nInput>(out, in, gauge, daggerFlag, Lx, Ly, Lz, Lt, grid_x, grid_y, grid_z, grid_t,
                                          parity, blockSize, this->dslashParam_.streams.stream);
}

template <typename _Float, int _nColor, int _nInput>
void WilsonDslash<_Float, _nColor, _nInput>::preApply() {}

template <typename _Float, int _nColor, int _nInput>
void WilsonDslash<_Float, _nColor, _nInput>::postApply() {}

template <typename _Float, int _nColor, int _nInput>
static void ApplyWilsonCore(DslashParam& dslashParam) {
  WilsonDslash<_Float, _nColor, _nInput> wilsonDslash(dslashParam);
  wilsonDslash.preApply();
  wilsonDslash.apply();
  wilsonDslash.postApply();
}

template <typename _Float, int _nColor>
inline void instantiateNumInput(const QcuLattDesc<Ndim>& lattDesc, const QcuProcDesc<Ndim>& procDesc, int parity,
                                int daggerFlag, double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* gauge,
                                int nInputs, ExecutionStreams& gpuStreams) {
  DslashParam dslashParam(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut, gauge,
                          gpuStreams);
  switch (nInputs) {
    case 1: {
      ApplyWilsonCore<_Float, _nColor, 1>(dslashParam);
    } break;
    case 2: {
      ApplyWilsonCore<_Float, _nColor, 2>(dslashParam);
    } break;
    case 3: {
      ApplyWilsonCore<_Float, _nColor, 3>(dslashParam);
    } break;
    case 4: {
      ApplyWilsonCore<_Float, _nColor, 4>(dslashParam);
    } break;
    case 5: {
      ApplyWilsonCore<_Float, _nColor, 5>(dslashParam);
    } break;

    default:
      errorQcu(
          "Now we only support 1-5 "
          "inputs.\n");
      break;
  }
}

template <typename _Float>
inline void instantiateColor(const QcuLattDesc<Ndim>& lattDesc, const QcuProcDesc<Ndim>& procDesc, int parity,
                             int daggerFlag, double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* pGauge,
                             int nColors, int nInputs, ExecutionStreams& gpuStreams) {
  switch (nColors) {
    case 3: {
      instantiateNumInput<_Float, 3>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut,
                                     pGauge, nInputs, gpuStreams);
    } break;
    case 4: {
      instantiateNumInput<_Float, 4>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut,
                                     pGauge, nInputs, gpuStreams);
    } break;
    case 5: {
      instantiateNumInput<_Float, 5>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut,
                                     pGauge, nInputs, gpuStreams);
    } break;
    case 6: {
      instantiateNumInput<_Float, 5>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut,
                                     pGauge, nInputs, gpuStreams);
    } break;
    case 7: {
      instantiateNumInput<_Float, 5>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut,
                                     pGauge, nInputs, gpuStreams);
    } break;
    case 8: {
      instantiateNumInput<_Float, 5>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut,
                                     pGauge, nInputs, gpuStreams);
    } break;

    default:
      errorQcu(
          "Now we only support 3-8 "
          "colors.\n");
      break;
  }
}

inline void instantiatePrec(const QcuLattDesc<Ndim>& lattDesc, const QcuProcDesc<Ndim>& procDesc, int parity,
                            int daggerFlag, double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* pGauge,
                            QCU_PRECISION floatPrecision, int nColors, int nInputs, ExecutionStreams& gpuStreams) {
  switch (floatPrecision) {
    // case QCU_HALF_PRECISION: {
    //   instantiateColor<half>(lattDesc,
    //   procDesc, parity, daggerFlag, kappa,
    //   mrhsFermionIn, mrhsFermionOut,
    //   pGauge, nColors, nInputs);
    // } break;
    case QCU_SINGLE_PRECISION: {
      instantiateColor<float>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut, pGauge,
                              nColors, nInputs, gpuStreams);
    } break;
    case QCU_DOUBLE_PRECISION: {
      instantiateColor<double>(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut, pGauge,
                               nColors, nInputs, gpuStreams);
    } break;

    default:
      errorQcu(
          "Now we only support single and "
          "double precision.\n");
      break;
  }
}
// void instantiateDslash(const QcuLattDesc<Nd>& lattDesc, const QcuProcDesc<Nd>& procDesc, int parity, int daggerFlag,
//                        double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* gauge, DSLASH_TYPE dtype,
//                        QCU_PRECISION floatPrecision, int nColors, int nInputs);
void instantiateDslash(const QcuLattDesc<Nd>& lattDesc, const QcuProcDesc<Nd>& procDesc, int parity, int daggerFlag,
                       double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* gauge, DSLASH_TYPE dtype,
                       QCU_PRECISION floatPrecision, int nColors, int nInputs, ExecutionStreams& gpuStreams) {
  switch (dtype) {
    case DSLASH_WILSON: {
      instantiatePrec(lattDesc, procDesc, parity, daggerFlag, kappa, mrhsFermionIn, mrhsFermionOut, gauge,
                      floatPrecision, nColors, nInputs, gpuStreams);
    } break;
    default:
      errorQcu(
          "Now we only support Wilson "
          "dslash.\n");
      break;
  }
}

}  // namespace qcu