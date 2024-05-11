#include "qcu_dirac.h"

namespace qcu {

// class DiracWilson : Dirac {
//  public:
void DiracWilson::Dslash(void* outputMRHS, void* inputMRHS, int parity, int daggerFlag) {
  // instantiateDslash(const QcuLattDesc<Nd>& lattDesc, const QcuProcDesc<Nd>& procDesc, int parity, int daggerFlag,
  //                      double kappa, void* mrhsFermionIn, void* mrhsFermionOut, void* gauge, DSLASH_TYPE dtype,
  //                      QCU_PRECISION floatPrecision, int nColors, int nInputs)
  DSLASH_TYPE dtype = DSLASH_WILSON;
  instantiateDslash(lattDesc_, procDesc_, parity, daggerFlag, kappa_, inputMRHS, outputMRHS, gauge_, dtype,
                    floatPrecision_, nColors_, nInputs_, streams_);
}

}  // namespace qcu