#ifndef KQED_SYMXY0_H
#define KQED_SYMXY0_H

// computes all kernels under full symmetrisation
__device__
void
compute_all_kernels_SYMXY0_v2( const double xv[4] ,
			       const double yv[4] ,
			       const struct QED_kernel_temps t ,
			       struct Kernels *K ) ;

// computes all kernels under full symmetrisation
__device__
void
compute_all_kernels_SYMXY0( const double xv[4] ,
			    const double yv[4] ,
			    const struct QED_kernel_temps t ,
			    struct QED_Kernels *K ) ;

#endif
