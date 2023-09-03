#ifndef KQED_H
#define KQED_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// alignment block
#define KQED_ALIGN

typedef enum { d0cb = 0 , d1cb = 1 , d2cb = 2 , d3cb = 3 } NDCB ;
typedef enum { QG0 , QG1 , QG2 , QG3 , QL4 , QL2 ,
	       dxQG0 , dxQG1 , dxQG2 , dxQG3 , dxQL4 , dxQL2 ,
	       d2xQG2 , d2xQG3  } FFidx ;
enum { YY = 0 , G0dy = 1 , G0dx = 2 , Gl2 = 3 , Gl2dy = 4 , Gl21 = 5 , Gl21dy = 6 , Gl3 = 7 , Gl3dy = 8 , G21 = 9 , G21dy = 10 , G22A = 11 ,  G22B = 12 , G22Ady = 13 , G22Bdy = 14 , G3A = 15 , G3B = 16 , G3Ady = 17 , G3Bdy = 18 , G31A = 19 , G31B = 20 , G31Ady = 21 , G31Bdy = 22 } ;
enum { alpha0dx_0p = 0 , alpha0_1p = 1 , alpha3_0p = 2 , beta2_1p = 3 , alpha3_1p = 4 , alpha1_0p = 5 , beta4_1p = 6 , alpha1_1p = 7 , alpha1dx_1p = 8 , alpha3dx_0p = 9 , alpha3dxdx_0p = 10 , beta2dx_1p = 11 , alpha3dx_1p = 12 , alpha1dx_0p = 13 } ;

struct Grid_coeffs {
  int Nffa ;
  float ****Ffp ;
  float ****Ffm ;
  double *XX ;
  double *YY ;
  double xstp ;
  double ystp ;
  int nstpx ;
  int nstpy ;
  double **TX ; // will be allocated to [23][NYTAY]
  int NtayX ;
  int NY_tay ;
  double **TY ; // will be allocated to [14][NXTAY]
  int NtayY ;
  int NX_tay ;
  int *nfx ;
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



#endif
