#ifndef KQED_SYMXY_H
#define KQED_SYMXY_H

// using our new data struct
__device__
void
compute_all_kernels_SYMXY_v2( const double xv[4] ,
			      const double yv[4] ,
			      const struct QED_kernel_temps t ,
			      struct Kernels *K ) ;

// computes all kernels symmetrising x and y
__device__
void
compute_all_kernels_SYMXY( const double xv[4] ,
			   const double yv[4] ,
			   const struct QED_kernel_temps t ,
			   struct QED_Kernels *K ) ;

#endif
