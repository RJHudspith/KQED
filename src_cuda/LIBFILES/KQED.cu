/**
 * Attempt at a single-TU ("unity") build.
 */

// for a unity build
#define KQED_PRIVATE static

#include "cheby.cu"
#include "getff-new.cu"
#include "chnr_dS.cu"
#include "chnr_dT.cu"
#include "chnr_dV.cu"
#include "Tabd.cu"

#include "QED_kernel.cu"
#include "QED_kernel_xy0.cu"
#include "kernels.cu"
#include "con_kernel.cu"
#include "sub_kernel.cu"
#include "all_kernels.cu"
#include "SYMXY.cu"
#include "SYMXY0.cu"

// #include "io.cu"
// #include "init.cu"
#include "pi_pert.cu"
