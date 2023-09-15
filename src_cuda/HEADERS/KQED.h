#ifndef KQED_H
#define KQED_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// is there a better way to include this!?
#include "config.h"

// alignment block
#define KQED_ALIGN

// with my version of gcc and -std=c99 this is not defined
#ifndef M_PI
  #define M_PI (3.14159265358979323846)
#endif

// if we have OMP we just use their routines
#if (defined HAVE_OMP_H) && (defined _OPENMP)
  #include <omp.h>
#else
  #define omp_get_thread_num()(0)
  #define omp_get_max_threads()(1)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// CUDA guards
#ifndef __device__
#define __device__
#endif
#ifndef __host__
#define __host__
#endif


// coupling constant of QED appears in pi_pert.c and amu_for_lattice
static const double AlfQED = 1.0/137.035999;

typedef enum { d0cb = 0 , d1cb = 1 , d2cb = 2 , d3cb = 3 } NDCB ;
typedef enum { QG0 = 0, QG1 , QG2 , QG3 , QL4 , QL2 ,
	       dxQG0 , dxQG1 , dxQG2 , dxQG3 , dxQL4 , dxQL2 ,
	       d2xQG2 , d2xQG3  } FFidx ;
enum { YY = 0 , G0dy = 1 , G0dx = 2 , Gl2 = 3 , Gl2dy = 4 , Gl21 = 5 , Gl21dy = 6 , Gl3 = 7 , Gl3dy = 8 , G21 = 9 , G21dy = 10 , G22A = 11 ,  G22B = 12 , G22Ady = 13 , G22Bdy = 14 , G3A = 15 , G3B = 16 , G3Ady = 17 , G3Bdy = 18 , G31A = 19 , G31B = 20 , G31Ady = 21 , G31Bdy = 22 } ;
enum { alpha0dx_0p = 0 , alpha0_1p = 1 , alpha3_0p = 2 , beta2_1p = 3 , alpha3_1p = 4 , alpha1_0p = 5 , beta4_1p = 6 , alpha1_1p = 7 , alpha1dx_1p = 8 , alpha3dx_0p = 9 , alpha3dxdx_0p = 10 , beta2dx_1p = 11 , alpha3dx_1p = 12 , alpha1dx_0p = 13 } ;

#define TX_LEN 23
#define TY_LEN 14
#define NYTAY 810
#define NXTAY 100

// NOTE: pointers are all expected to be device ptrs
struct Grid_coeffs {
  int Nffa ;
  float *Ffp ; // 4D array [Nffa][nstpx][nstpy][nfx_max]
  float *Ffm ; // 4D array [Nffa][nstpx][nstpy][nfx_max]
  double *XX ; // 1D array [nstpx]
  double *YY ; // 1D array [nstpy]
  double xstp ;
  double ystp ;
  int nstpx ;
  int nstpy ;
  double *TX ; // 2D array [NtayY][NY_tay]
  int NtayY ; // = TX_LEN
  int NY_tay ; // = NYTAY
  double *TY ; // 2D array [NtayX][NX_tay]
  int NtayX ; // = TY_LEN
  int NX_tay ; // = NXTAY
  int nfx_max;
  int *nfx ; // 1D array [nstpx]
} ;

__device__ inline static double* getTX(const struct Grid_coeffs* Grid, int i) {
  return &Grid -> TX[i * Grid -> NY_tay];
}
__device__ inline static double* getTY(const struct Grid_coeffs* Grid, int i) {
  return &Grid -> TY[i * Grid -> NX_tay];
}
__device__ inline static float* getFfp(const struct Grid_coeffs* Grid, int i, int j, int k) {
  return &Grid -> Ffp[((i * Grid -> nstpx + j) * Grid -> nstpy + k) * Grid -> nfx_max];
}
__device__ inline static float* getFfm(const struct Grid_coeffs* Grid, int i, int j, int k) {
  return &Grid -> Ffm[((i * Grid -> nstpx + j) * Grid -> nstpy + k) * Grid -> nfx_max];
}

// Allocations occur in this
// NOTE: pointers are all expected to be device ptrs
struct QED_kernel_temps {
  struct Grid_coeffs Grid ;
  int *G8 ;
} ;


struct STV {
  double Sxv[4] KQED_ALIGN ;
  double Syv[4] KQED_ALIGN ;
  double Txv[4][4][4] KQED_ALIGN ;
  double Tyv[4][4][4] KQED_ALIGN ;
  double Vv[4][4][4] KQED_ALIGN ;
} ;
struct Kernel {
  double xy[6][4][4][4] KQED_ALIGN ;
  double yx[6][4][4][4] KQED_ALIGN ;
} ;
struct Kernels {
  double L0[6][4][4][4] KQED_ALIGN ;
  double L1[6][4][4][4] KQED_ALIGN ;
  double L2[6][4][4][4] KQED_ALIGN ;
  double L3[6][4][4][4] KQED_ALIGN ;
} ;
struct QED_Kernels {
  struct Kernel L0 ;
  struct Kernel L1 ;
  struct Kernel L2 ;
  struct Kernel L3 ;
} ;
struct QED_Con_Kernels {
  struct Kernels A ; // L^i(x,y)
  struct Kernels B ; // L^i(x,x-y)
  struct Kernels C ; // combination L^i(x,y)+L^i_{mu<>nu}(y,x)-L^i(x,x-y)
} ;

struct intprecomp {
  double A , B , C1 , C2 , D , lA ;
  size_t idx ;
} ;
struct invariants {
  double x,cb,y,xsq, ysq, xdoty;
  double cborig , yorig ;
  double xmysq ;
  double xmy ;
  bool flag2 ;
  struct intprecomp INVx , INVy ;
} ;


// typical defines for the main example if we need more functionality
// just add more headers here
#include "all_kernels.h"
#include "con_kernel.h"
#include "sub_kernel.h"
#include "GLU_timer.h"
#include "init.h"
#include "kernels.h"
#include "pi_pert.h"
#include "SYMXY.h"
#include "SYMXY0.h"

#ifdef __cplusplus
} // extern "C"
#endif


#endif
