#pragma once

#include "qcu.h"
#include "qcuDesc/mrhs_desc.h"
#include "qcuDesc/qcu_desc.h"
#include "qcu_enum.h"
#include "qcu_kernel_param.cuh"

namespace qcu {

static constexpr int Ndim = 4;

struct DslashParam {
  const QcuProcDesc<Ndim>& procDesc;
  const QcuLattDesc<Ndim>& lattDesc;
  int parity;
  int daggerFlag;
  double kappa;

  void* fermionIn_MRHS;
  void* fermionOut_MRHS;

  void* gauge;
  ExecutionStreams& streams;
  DslashParam(const QcuLattDesc<Ndim>& pLattDesc, const QcuProcDesc<Ndim>& pProcDesc, int pParity, int pDaggerFlag,
              double pKappa, void* pFermionIn, void* pFermionOut, void* pGauge, ExecutionStreams& pStreams)
      : lattDesc(pLattDesc),
        procDesc(pProcDesc),
        parity(pParity),
        daggerFlag(pDaggerFlag),
        kappa(pKappa),
        fermionIn_MRHS(pFermionIn),
        fermionOut_MRHS(pFermionOut),
        gauge(pGauge),
        streams(pStreams) {}
};

template <typename _Float, int _nColor, int _nInput>
class Dslash {
 protected:
  const DslashParam& dslashParam_;

 public:
  Dslash(const DslashParam& dslashParam) : dslashParam_(dslashParam) {}
  virtual ~Dslash() {}
  virtual void apply() = 0;
  virtual void preApply() = 0;
  virtual void postApply() = 0;
};

// template <template <int, int> class _DParam, typename _Float, int _nColor, int _nInput>
template <typename _Float, int _nColor, int _nInput>
class WilsonDslash : public Dslash<_Float, _nColor, _nInput> {
 public:
  WilsonDslash(const DslashParam& dslashParam) : Dslash<_Float, _nColor, _nInput>(dslashParam) {}
  virtual ~WilsonDslash() {}
  virtual void apply();
  virtual void preApply();
  virtual void postApply();
};

void instantiateDslash(const QcuLattDesc<Nd>& lattDesc, const QcuProcDesc<Nd>& procDesc, int parity, int daggerFlag,
                       double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* gauge, DSLASH_TYPE dtype,
                       QCU_PRECISION floatPrecision, int nColors, int nInputs, ExecutionStreams& streams);
// template
}  // namespace qcu