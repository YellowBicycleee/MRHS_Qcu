#pragma once

#include <cassert>

#include "qcu.h"

// in every process, lattice size desc
template <int _Ndim>
struct QcuLattDesc {
  int lattSize[_Ndim];

  QcuLattDesc(int x, int y, int z, int t) {
    lattSize[X_DIM] = x;
    if (Y_DIM < _Ndim) lattSize[Y_DIM] = y;
    if (Z_DIM < _Ndim) lattSize[Z_DIM] = z;
    if (T_DIM < _Ndim) lattSize[T_DIM] = t;
  }
  QcuLattDesc(QcuParam *param) {
    lattSize[X_DIM] = param->lattice_size[X_DIM];
    if (Y_DIM < _Ndim) lattSize[Y_DIM] = param->lattice_size[Y_DIM];
    if (Z_DIM < _Ndim) lattSize[Z_DIM] = param->lattice_size[Z_DIM];
    if (T_DIM < _Ndim) lattSize[T_DIM] = param->lattice_size[T_DIM];
  }

  int X() { return lattSize[X_DIM]; }
  int Y() {
    if (_Ndim > Y_DIM)
      return lattSize[Y_DIM];
    else {
      printf("file %s line %d dim = %d\n", __FILE__, __LINE__, _Ndim);
      assert(0);
    }
  }
  int Z() {
    if (_Ndim > Z_DIM)
      return lattSize[Z_DIM];
    else {
      printf("file %s line %d dim = %d\n", __FILE__, __LINE__, _Ndim);
      assert(0);
    }
  }
  int T() {
    if (_Ndim > T_DIM)
      return lattSize[T_DIM];
    else {
      printf("file %s line %d dim = %d\n", __FILE__, __LINE__, _Ndim);
      assert(0);
    }
  }
};

template <int _Ndim>
struct QcuProcDesc {
  int procSize[_Ndim];

  QcuProcDesc(int x, int y, int z, int t) {
    procSize[X_DIM] = x;
    if (Y_DIM < _Ndim) procSize[Y_DIM] = y;
    if (Z_DIM < _Ndim) procSize[Z_DIM] = z;
    if (T_DIM < _Ndim) procSize[T_DIM] = t;
  }
  QcuProcDesc(QcuGrid *grid) {
    procSize[X_DIM] = grid->grid_size[X_DIM];
    if (Y_DIM < _Ndim) procSize[Y_DIM] = grid->grid_size[Y_DIM];
    if (Z_DIM < _Ndim) procSize[Z_DIM] = grid->grid_size[Z_DIM];
    if (T_DIM < _Ndim) procSize[T_DIM] = grid->grid_size[T_DIM];
  }

  int X() { return procSize[X_DIM]; }
  int Y() {
    if (_Ndim > Y_DIM)
      return procSize[Y_DIM];
    else {
      printf("file %s line %d dim = %d\n", __FILE__, __LINE__, _Ndim);
      assert(0);
    }
  }

  int Z() {
    if (_Ndim > Z_DIM)
      return procSize[Z_DIM];
    else {
      printf("file %s line %d dim = %d\n", __FILE__, __LINE__, _Ndim);
      assert(0);
    }
  }

  int T() {
    if (_Ndim > T_DIM)
      return procSize[T_DIM];
    else {
      printf("file %s line %d dim = %d\n", __FILE__, __LINE__, _Ndim);
      assert(0);
    }
  }
};
