#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QcuParam_s {
  int lattice_size[4];
} QcuParam;

typedef struct QcuGrid_t {
  int grid_size[4];
} QcuGrid_t;

typedef QcuGrid_t QcuGrid;

// void initGridSize(QcuGrid_t *grid, QcuParam *p_param, void *gauge, void *fermion_in, void *fermion_out);
// void initGridSize(QcuGrid_t *grid, QcuParam *p_param);
void initGridSize(QcuGrid_t *grid, QcuParam *p_param, int nColors);
void pushBackFermions(void *fermionOut, void *fermionIn);  // 押入输入输出向量
void startDslash(int parity, int daggerFlag);
void loadQcuGauge(void *gauge);
void getDslash(int dslashType, void *gauge, double mass, int nColors, int nInputs, int floatPrecision);
// void loadQcuGauge(void *gauge, QcuParam *param); // 废弃

void dslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int parity);           // 废弃
void fullDslashQcu(void *fermion_out, void *fermion_in, void *gauge, QcuParam *param, int dagger_flag);  // 废弃
// TODO
void cg_inverter(void *x_vector, void *b_vector, void *gauge, QcuParam *param, double p_max_prec, double p_kappa);

#ifdef __cplusplus
}
#endif
