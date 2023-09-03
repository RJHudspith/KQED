#ifndef QED_KERNEL_H
#define QED_KERNEL_H

__device__
int
kernelQED( const double xv[4] ,
	   const double yv[4] ,
	   struct QED_kernel_temps t ,
	   double kerv[6][4][4][4] ) ;

#endif
