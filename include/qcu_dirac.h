#pragma once
#include "qcd/qcu_dslash.h"
#include "qcuDesc/qcu_desc.h"
#include "qcu_kernel_param.cuh"
namespace qcu {

class Dirac {
 protected:
  void* gauge_;
  double kappa_;
  double mass_;
  int nColors_;
  //   int nInputs_;
  QCU_PRECISION floatPrecision_;
  const QcuLattDesc<Nd>& lattDesc_;
  const QcuProcDesc<Nd>& procDesc_;
  DSLASH_TYPE dslashType_;

  ExecutionStreams& streams_;

 public:
  Dirac(void* gauge, double mass, int nColors, QCU_PRECISION floatPrecision, const QcuLattDesc<Nd>& lattDesc,
        const QcuProcDesc<Nd>& procDesc, DSLASH_TYPE dslashType, ExecutionStreams& streams)
      : gauge_(gauge),
        mass_(mass),
        kappa_(1.0 / (2.0 * (4.0 + mass))),
        nColors_(nColors),
        // nInputs_(nInputs),
        floatPrecision_(floatPrecision),
        lattDesc_(lattDesc),
        procDesc_(procDesc),
        dslashType_(dslashType),
        streams_(streams) {}

  virtual ~Dirac() {}
  virtual void Dslash(void* outputMRHS, void* inputMRHS, int parity, int daggerFlag, int nInputs) = 0;
  virtual void DslashXpay() = 0;
};

class DiracWilson : public Dirac {
 public:
  DiracWilson(void* gauge, double mass, int nColors, QCU_PRECISION floatPrecision, const QcuLattDesc<Nd>& lattDesc,
              const QcuProcDesc<Nd>& procDesc, DSLASH_TYPE dslashType, ExecutionStreams& streams)
      : Dirac(gauge, mass, nColors, floatPrecision, lattDesc, procDesc, dslashType, streams) {}
  ~DiracWilson() {}
  // TODO
  virtual void Dslash(void* outputMRHS, void* inputMRHS, int parity, int daggerFlag, int nInputs);
  // TODO
  virtual void DslashXpay() {}
};

}  // namespace qcu