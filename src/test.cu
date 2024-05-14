#include <iostream>
#include "qcu.h"
#include "qcu_macro.cuh"
#include "qcu_enum.h"
#include "storage/qcu_storage.h"
using namespace std;


int main () {
    void* fermionIn;
    void* fermionOut;
    void* gauge;
    int Lt = 8;
    int Lz = 8;
    int Ly = 8;
    int Lx = 8;
    int vol = Lx * Ly * Lz * Lt;
    int Nc = 3;
    printf("****colorSpinor len = %d, gauge %d********\n", 
	2 * sizeof(double) * vol / 2 * Ns * Nc,
	2 * sizeof(double) * vol * 4 * Nc * Nc
    );
    CHECK_CUDA(cudaMalloc(&fermionIn, 2 * sizeof(double) * vol / 2 * Ns * Nc));
    CHECK_CUDA(cudaMalloc(&fermionOut, 2 * sizeof(double) * vol / 2 * Ns * Nc));
    CHECK_CUDA(cudaMalloc(&gauge, 2 * sizeof(double) * vol * Nc * Nc * Nd));

    // void initGridSize(QcuGrid_t *grid, QcuParam *p_param, int nColors);
    QcuGrid_t grid;// = {.grid_size[0] = 1, .grid_size[1] = 1, .grid_size[2] = 1, .grid_size[3] = 1 };
    grid.grid_size[0] = 1;
    grid.grid_size[1] = 1;
    grid.grid_size[2] = 1;
    grid.grid_size[3] = 1;
    QcuParam param; // = {.lattice_size[0] = Lx, .lattice_size[1] = Ly, .lattice_size[2] = Lz, .lattice_size[3] = Lt };
    param.lattice_size[0] = Lx;
    param.lattice_size[1] = Ly;
    param.lattice_size[2] = Lx;
    param.lattice_size[3] = Lt;
    for (int i = 0; i< 1; i++) {
        initGridSize(&grid, &param, 3);
        loadQcuGauge(gauge);
        int dslashType = 0;
        printf("===================DEBUG BEGIN=========\n");
        qcu::storage::shiftVectorStorage(fermionOut, fermionIn, TO_COALESCE, Lx / 2, Ly, Lz, Lt);

        printf("===================DEBUG END=========\n");
        getDslash(dslashType, -3.5, Nc, 1, 2);
        pushBackFermions(fermionOut, fermionIn);
        startDslash(0, 0);
    }
    CHECK_CUDA(cudaFree(fermionIn));
    CHECK_CUDA(cudaFree(fermionOut));
    CHECK_CUDA(cudaFree(gauge));
    return 0;
}

