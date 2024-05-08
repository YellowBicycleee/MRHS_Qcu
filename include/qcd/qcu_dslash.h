#pragma once

#include "qcu.h"
#include "qcuDesc/qcuDesc.h"
#include "qcu_enum.h"

struct DslashParam {
  const QcuProcDesc<4>& procDesc;
  const QcuLattDesc<4>& lattDesc;
  int parity;
  DslashParam(const QcuLattDesc& pLattDesc, const QcuProcDesc& pProcessDesc, int pParity)
      : lattDesc(pLattDesc), processDesc(pProcessDesc), parity(pParity) {}
};


template <int _nColor, QcuPrecondition _precondition = QCU_EVEN_ODD_PRECONDITION>
class Dslash {
  const DslashParam& dslashParam_;
 public:
  Dslash(const DslashParam& dslashParam) : dslashParam_(dslashParam) {}

  virtual void apply();
  virtual void preApply() {}
  virtual void postApply() {}
}