#ifndef KQED_CON_KERNEL_H
#define KQED_CON_KERNEL_H

// average two con kernels
__device__
void
average_con_kernels( struct Kernels *K ,
		     const struct Kernels *W ) ;

// compute the kernel combinations needed for the connected contribution
__device__
void
compute_con_kernels( const double xv[4] ,
		     const double yv[4] ,
		     const struct QED_kernel_temps t ,
		     struct QED_Kernels *K ) ;

__device__
void
compute_con_kernels_v2( const double xv[4] ,
			const double yv[4] ,
			const struct QED_kernel_temps t ,
			struct Kernels *K ) ;

__device__
void
compute_con_kernelsM_L2( const double M[4] ,
			 const double xv[4] ,
			 const double yv[4] ,
			 const struct QED_kernel_temps t ,
			 struct Kernels *K ) ;

#endif
