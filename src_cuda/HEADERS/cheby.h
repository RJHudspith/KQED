#ifndef CHEBY_H
#define CHEBY_H

// evaluation of cheby polynomial
__device__ KQED_PRIVATE
double
chebUsum( const int nk ,
	  const double x ,
	  const double *co ) ;

// first order deriv of cheby polynomial
__device__ KQED_PRIVATE
double
dchebUsum( const int nk,
	   const double x,
	   const double *co) ;

// second order deriv of cheby polynomial
__device__ KQED_PRIVATE
double
ddchebUsum( const int nk,
	    const double x,
	    const double *co) ;

// third-order deriv of cheby polynomial
__device__ KQED_PRIVATE
double
dddchebUsum( const int nk ,
	     const double x ,
	     const double *co ) ;

#endif
