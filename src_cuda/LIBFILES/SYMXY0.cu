/**
   @file SYMXY0.c
   @brief performs the full symmetrisation of the kernel

   It is full because it covers all the permutations of indices

   (Dropping \rho \sigma indices...) 
   
   L^{SYMXY0}_{\mu\nu\lambda}(x,y) = 
      ( +L_{\mu\nu\lambda}(x,y) + L_{\nu\mu\lambda}(y,x) 
        -L_{\lambda\nu\mu}(x,x-y) - L_{\nu\lambda\nu}(x-y,x)
        -L_{\lambda\mu\nu}(y,y-x) - L_{\mu\lambda\nu}(y-x,y) )/6

   And sets

   L^{SYMXY0}_{\mu\nu\lambda}(x,y) = L^{SYMXY0}_{\nu\mu\lambda}(y,x) = 
 */
#include "KQED.h"

#include "kernels.h"


// computes all kernels under full symmetrisation
__device__
void
compute_all_kernels_SYMXY0_v2( const double xv[4] ,
			       const double yv[4] ,
			       const struct QED_kernel_temps t ,
			       struct Kernels *K )
{
  // compute x-y
  const double xmyv[4] = { xv[0]-yv[0] , xv[1]-yv[1] ,
			   xv[2]-yv[2] , xv[3]-yv[3] } ;
  // need the kernels at L(x-y,x), L(y-x,y), L(x,x-y), and L(y,y-x)
  double Lxmy_x[6][4][4][4] , Ly_x[6][4][4][4] ;
  QED_kernel_L0( xv   , yv , t , K -> L0 ) ;
  QED_kernel_L0( yv   , xv , t , Ly_x ) ;
  QED_kernel_L0( xmyv , xv , t , Lxmy_x ) ;
  // do the y-x computation
  const double ymxv[4] = { -xmyv[0] , -xmyv[1] , -xmyv[2] , -xmyv[3] } ;
  double Lymx_y[6][4][4][4] , Lx_xmy[6][4][4][4] , Ly_ymx[6][4][4][4] ;
  QED_kernel_L0( ymxv , yv   , t , Lymx_y ) ;
  QED_kernel_L0( xv   , xmyv , t , Lx_xmy ) ;
  QED_kernel_L0( yv   , ymxv , t , Ly_ymx ) ;
  
  // and these temporaries: L(x,0) L(y,0) L(0,x) L(0,y)
  // and also L(x-y,0) and L(0,x-y)
  // we can rebuild L(y-x,0) from -L(x-y) and likewise L(0,y-x)
  // so we don't compute them
  double Lx_0[6][4][4][4]   , L0_x[6][4][4][4]   ;
  double Ly_0[6][4][4][4]   , L0_y[6][4][4][4]   ;
  double Lxmy_0[6][4][4][4] , L0_xmy[6][4][4][4] ;  
  const double zero[4] = { 0. , 0. , 0. , 0. } ;
  QED_kernel_L0( xv   , zero , t , Lx_0 )   ;
  QED_kernel_L0( zero , xv   , t , L0_x )   ;
  QED_kernel_L0( yv   , zero , t , Ly_0 )   ;
  QED_kernel_L0( zero , yv   , t , L0_y )   ;
  QED_kernel_L0( xmyv , zero , t , Lxmy_0 ) ;
  QED_kernel_L0( zero , xmyv , t , L0_xmy ) ;

  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // symmetrised L0
	  K -> L0[rhosig][mu][nu][lambda] =
	    ( +K -> L0[rhosig][mu][nu][lambda] 
	      +Ly_x[rhosig][nu][mu][lambda]
	      -Lxmy_x[rhosig][nu][lambda][mu]
	      -Lymx_y[rhosig][mu][lambda][nu]
	      -Lx_xmy[rhosig][lambda][nu][mu]
	      -Ly_ymx[rhosig][lambda][mu][nu] )/6 ;
	  
	  // this is a common temporary in L1 and L2
	  // (also in L3 but that has a substantial contribution from L2)
	  register const double ApC =
	    +Lx_0[rhosig][mu][nu][lambda] - Lx_0[rhosig][lambda][nu][mu]
	    +Ly_0[rhosig][nu][mu][lambda] - Ly_0[rhosig][lambda][mu][nu] ;
	      
	  // symmetrised L1
	  K -> L1[rhosig][mu][nu][lambda] =
	    K -> L0[rhosig][mu][nu][lambda] +
	    ( -ApC
	      +Lx_0[rhosig][lambda][mu][nu]    - Lx_0[rhosig][mu][lambda][nu]
	      +Ly_0[rhosig][lambda][nu][mu]    - Ly_0[rhosig][nu][lambda][mu]
	      +Lxmy_0[rhosig][nu][mu][lambda]  - Lxmy_0[rhosig][mu][nu][lambda]
	      +Lxmy_0[rhosig][nu][lambda][mu]  - Lxmy_0[rhosig][mu][lambda][nu] )/12 ;
	  
	  // symmetrised L2
	  K -> L2[rhosig][mu][nu][lambda] =
	    K -> L0[rhosig][mu][nu][lambda] -
	    ( +ApC
	      +L0_x[rhosig][nu][mu][lambda]    - L0_x[rhosig][nu][lambda][mu]
	      +L0_y[rhosig][mu][nu][lambda]    - L0_y[rhosig][mu][lambda][nu]
	      +Lxmy_0[rhosig][mu][lambda][nu]  - Lxmy_0[rhosig][nu][lambda][mu]  
	      +L0_xmy[rhosig][lambda][mu][nu]  - L0_xmy[rhosig][lambda][nu][mu] )/6. ;

	  // symmetrised L3 is basically L2 with some other junk
	  K -> L3[rhosig][mu][nu][lambda] =
	    K -> L2[rhosig][mu][nu][lambda] +
	    ( +L0_x[rhosig][mu][nu][lambda]    - L0_x[rhosig][lambda][nu][mu]
	      +L0_y[rhosig][nu][mu][lambda]    - L0_y[rhosig][lambda][mu][nu]
	      +L0_xmy[rhosig][mu][lambda][nu]  - L0_xmy[rhosig][nu][lambda][mu] )/6. ;	  
	}
      }
    }
  }
  
  return ;
}

// computes all kernels under full symmetrisation
__device__
void
compute_all_kernels_SYMXY0( const double xv[4] ,
			    const double yv[4] ,
			    const struct QED_kernel_temps t ,
			    struct QED_Kernels *K )
{
  // compute x-y
  const double xmyv[4] = { xv[0]-yv[0] , xv[1]-yv[1] ,
			   xv[2]-yv[2] , xv[3]-yv[3] } ;
  // need the kernels at L(x-y,x), L(y-x,y), L(x,x-y), and L(y,y-x)
  double Lxmy_x[6][4][4][4] ;
  QED_kernel_L0( xv   , yv , t , K -> L0.xy ) ;
  QED_kernel_L0( yv   , xv , t , K -> L0.yx ) ;
  QED_kernel_L0( xmyv , xv , t , Lxmy_x ) ;
  // do the y-x computation
  const double ymxv[4] = { -xmyv[0] , -xmyv[1] , -xmyv[2] , -xmyv[3] } ;
  double Lymx_y[6][4][4][4] , Lx_xmy[6][4][4][4] , Ly_ymx[6][4][4][4] ;
  QED_kernel_L0( ymxv , yv   , t , Lymx_y ) ;
  QED_kernel_L0( xv   , xmyv , t , Lx_xmy ) ;
  QED_kernel_L0( yv   , ymxv , t , Ly_ymx ) ;
  
  // and these temporaries: L(x,0) L(y,0) L(0,x) L(0,y)
  // and also L(x-y,0) and L(0,x-y)
  // we can rebuild L(y-x,0) from -L(x-y) and likewise L(0,y-x)
  // so we don't compute them
  double Lx_0[6][4][4][4]   , L0_x[6][4][4][4]   ;
  double Ly_0[6][4][4][4]   , L0_y[6][4][4][4]   ;
  double Lxmy_0[6][4][4][4] , L0_xmy[6][4][4][4] ;  
  const double zero[4] = { 0. , 0. , 0. , 0. } ;
  QED_kernel_L0( xv   , zero , t , Lx_0 )   ;
  QED_kernel_L0( zero , xv   , t , L0_x )   ;
  QED_kernel_L0( yv   , zero , t , Ly_0 )   ;
  QED_kernel_L0( zero , yv   , t , L0_y )   ;
  QED_kernel_L0( xmyv , zero , t , Lxmy_0 ) ;
  QED_kernel_L0( zero , xmyv , t , L0_xmy ) ;


  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // symmetrised L0
	  K -> L0.xy[rhosig][mu][nu][lambda] =
	    ( +K -> L0.xy[rhosig][mu][nu][lambda] 
	      +K -> L0.yx[rhosig][nu][mu][lambda]
	      -Lxmy_x[rhosig][nu][lambda][mu]
	      -Lymx_y[rhosig][mu][lambda][nu]
	      -Lx_xmy[rhosig][lambda][nu][mu]
	      -Ly_ymx[rhosig][lambda][mu][nu] )/6 ;
	  
	  // this is a common temporary in L1 and L2
	  // (also in L3 but that has a substantial contribution from L2)
	  register const double ApC =
	    +Lx_0[rhosig][mu][nu][lambda] - Lx_0[rhosig][lambda][nu][mu]
	    +Ly_0[rhosig][nu][mu][lambda] - Ly_0[rhosig][lambda][mu][nu] ;
	      
	  // symmetrised L1
	  K -> L1.xy[rhosig][mu][nu][lambda] =
	    K -> L0.xy[rhosig][mu][nu][lambda] +
	    ( -ApC
	      +Lx_0[rhosig][lambda][mu][nu]    - Lx_0[rhosig][mu][lambda][nu]
	      +Ly_0[rhosig][lambda][nu][mu]    - Ly_0[rhosig][nu][lambda][mu]
	      +Lxmy_0[rhosig][nu][mu][lambda]  - Lxmy_0[rhosig][mu][nu][lambda]
	      +Lxmy_0[rhosig][nu][lambda][mu]  - Lxmy_0[rhosig][mu][lambda][nu] )/12 ;
	  
	  // symmetrised L2
	  K -> L2.xy[rhosig][mu][nu][lambda] =
	    K -> L0.xy[rhosig][mu][nu][lambda] -
	    ( +ApC
	      +L0_x[rhosig][nu][mu][lambda]    - L0_x[rhosig][nu][lambda][mu]
	      +L0_y[rhosig][mu][nu][lambda]    - L0_y[rhosig][mu][lambda][nu]
	      +Lxmy_0[rhosig][mu][lambda][nu]  - Lxmy_0[rhosig][nu][lambda][mu]  
	      +L0_xmy[rhosig][lambda][mu][nu]  - L0_xmy[rhosig][lambda][nu][mu] )/6. ;

	  // symmetrised L3 is basically L2 with some other junk
	  K -> L3.xy[rhosig][mu][nu][lambda] =
	    K -> L2.xy[rhosig][mu][nu][lambda] +
	    ( +L0_x[rhosig][mu][nu][lambda]    - L0_x[rhosig][lambda][nu][mu]
	      +L0_y[rhosig][nu][mu][lambda]    - L0_y[rhosig][lambda][mu][nu]
	      +L0_xmy[rhosig][mu][lambda][nu]  - L0_xmy[rhosig][nu][lambda][mu] )/6. ;	  
	}
      }
    }
  }

  // and set Lyx to Lxy with a swap of mu and nu indices
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  K->L0.yx[rhosig][mu][nu][lambda] = K->L0.xy[rhosig][nu][mu][lambda] ;
	  K->L1.yx[rhosig][mu][nu][lambda] = K->L1.xy[rhosig][nu][mu][lambda] ;
	  K->L2.yx[rhosig][mu][nu][lambda] = K->L2.xy[rhosig][nu][mu][lambda] ;
	  K->L3.yx[rhosig][mu][nu][lambda] = K->L3.xy[rhosig][nu][mu][lambda] ;
	}
      }
    }
  }
  
  return ;
}
