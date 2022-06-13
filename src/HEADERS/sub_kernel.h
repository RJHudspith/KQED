#ifndef KQED_SUB_KERNEL_H
#define KQED_SUB_KERNEL_H

void
compute_sub_kernelsM_L2( const double M[4] ,
			 const double xv[4] ,
			 const double yv[4] ,
			 const struct QED_kernel_temps t ,
			 struct Kernels *K ) ;

#endif
