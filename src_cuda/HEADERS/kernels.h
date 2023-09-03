#ifndef KERNELS_H
#define KERNELS_H

__device__
void
QED_kernel_L0( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

__device__
void
QED_kernel_L1( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

__device__
void
QED_kernel_L2( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

__device__
void
QED_Mkernel_L2( const double M ,
		const double xv[4] ,
		const double yv[4] ,
		const struct QED_kernel_temps t ,
		double kerv[6][4][4][4] ) ;

__device__
void
QED_kernel_L3( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

#endif
