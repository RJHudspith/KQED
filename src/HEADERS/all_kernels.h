#ifndef KQED_ALL_KERNELS_H
#define KQED_ALL_KERNELS_H

void
average_all_kernels( struct Kernels *K ,
		     const struct Kernels W ) ;

void
average_all_QED_kernels( struct QED_Kernels *K ,
			 const struct QED_Kernels W ) ;

void
compute_all_kernels( const double xv[4] ,
		     const double yv[4] ,
		     const struct QED_kernel_temps t ,
		     struct QED_Kernels *K ) ;

void
compute_all_Mkernels( const double M[4] ,
		      const double xv[4] ,
		      const double yv[4] ,
		      const struct QED_kernel_temps t ,
		      struct QED_Kernels *K ) ;

void
compute_all_Mkernels_v2( const double M[4] ,
			 const double xv[4] ,
			 const double yv[4] ,
			 const struct QED_kernel_temps t ,
			 struct QED_Kernels *K ) ;

void
swap_munu_Lyx( struct QED_Kernels *K ) ;

#endif
