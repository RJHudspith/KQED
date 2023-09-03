/**
   @file SYMXY.c
   @brief symmetrises x and y arguments of the kernel

   L^{SYMXY}_{\rho\sigma;\mu\nu\lambda}(x,y)
      = ( L^{SYMXY}_{\rho\sigma;\mu\nu\lambda}(x,y)
            + L^{SYMXY}_{\rho\sigma;\nu\mu\lambda}(y,x) )/2
 */
#include "KQED.h"

#include "all_kernels.h"

// computes all kernels symmetrising x and y
__device__
void
compute_all_kernels_SYMXY_v2( const double xv[4] ,
			      const double yv[4] ,
			      const struct QED_kernel_temps t ,
			      struct Kernels *K )
{
  struct QED_Kernels W ;
  compute_all_kernels( xv , yv , t , &W ) ;
  
  // perform the symmetrisation L(x,y) = ( L(x,y) + L_{mu<->nu}(y,x))/2
  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  K->L0[rhosig][mu][nu][lambda] = ( W.L0.xy[rhosig][mu][nu][lambda] +
					    W.L0.yx[rhosig][nu][mu][lambda] )/2 ; 
	  K->L1[rhosig][mu][nu][lambda] = ( W.L1.xy[rhosig][mu][nu][lambda] +
					    W.L1.yx[rhosig][nu][mu][lambda] )/2 ;
	  K->L2[rhosig][mu][nu][lambda] = ( W.L2.xy[rhosig][mu][nu][lambda] +
					    W.L2.yx[rhosig][nu][mu][lambda] )/2 ;
	  K->L3[rhosig][mu][nu][lambda] = ( W.L3.xy[rhosig][mu][nu][lambda] +
					    W.L3.yx[rhosig][nu][mu][lambda] )/2 ; 
	}
      }
    }
  }
  return ;
}

// computes all kernels symmetrising x and y
__device__
void
compute_all_kernels_SYMXY( const double xv[4] ,
			   const double yv[4] ,
			   const struct QED_kernel_temps t ,
			   struct QED_Kernels *K )
{
  compute_all_kernels( xv , yv , t , K ) ;

  // perform the symmetrisation L(x,y) = ( L(x,y) + L_{mu<->nu}(y,x))/2
  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  K->L0.xy[rhosig][mu][nu][lambda] = ( +K->L0.xy[rhosig][mu][nu][lambda]
					       +K->L0.yx[rhosig][nu][mu][lambda] )/2 ; 
	  K->L1.xy[rhosig][mu][nu][lambda] = ( +K->L1.xy[rhosig][mu][nu][lambda]
					       +K->L1.yx[rhosig][nu][mu][lambda] )/2 ;
	  K->L2.xy[rhosig][mu][nu][lambda] = ( +K->L2.xy[rhosig][mu][nu][lambda]
					       +K->L2.yx[rhosig][nu][mu][lambda] )/2 ;
	  K->L3.xy[rhosig][mu][nu][lambda] = ( +K->L3.xy[rhosig][mu][nu][lambda]
					       +K->L3.yx[rhosig][nu][mu][lambda] )/2 ; 
	}
      }
    }
  }

    // and set Lyx to Lxy :: note that the indices are the same as we have symmetrised
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
