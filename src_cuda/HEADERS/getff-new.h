#ifndef GETFF_NEW_H
#define GETFF_NEW_H

#include "KQED.h"

// accesses V apparently
__device__ KQED_PRIVATE
double
accessv( const bool flag_hy , const bool use_y_derivs,
	 const int ix, const int iy,
	 const FFidx nm, const bool ndy , const NDCB ndcb ,
	 const double cb , const double y ,
	 const struct Grid_coeffs Grid ) ;

// performs search on a sorted list returning the closest coordinate
__device__ KQED_PRIVATE
int
find_ind( const double *arr , const double target ,
       const int lo , const int hi ) ;

__device__ KQED_PRIVATE
double
extractff( const FFidx nm, const bool ndy, const NDCB ndcb ,
	   const struct invariants Inv ,
	   const struct Grid_coeffs Grid ) ;

__device__ KQED_PRIVATE
void
extractff2( const FFidx nm, const NDCB ndcb ,
	    const struct invariants Inv ,
	    const struct Grid_coeffs Grid ,
	    double f[4] ) ;

// little helper function
__device__
static inline double
lerp( const double a ,
      const double T1 ,
      const double T2 )
{
  return a*T1 + (1.0-a)*T2 ;
}

// precompute all this business for x and y
__device__ KQED_PRIVATE
void
precompute_INV( struct intprecomp *INV ,
		 const double y ,
		 const double y1 ,
		 const double y2 ,
		 const size_t idx ) ;
__device__ KQED_PRIVATE
void
precompute_INVx( struct intprecomp *INVx ,
		 const double y ,
		 const double y1 ,
		 const double y2 ,
		 const size_t idx ) ;

// scalar product of two four-vectors
__device__
static inline double
SCALPROD( const double xv[4] ,
	  const double yv[4] )
{
  return xv[0]*yv[0] + xv[1]*yv[1] + xv[2]*yv[2] + xv[3]*yv[3] ;
}

// initialises the invariants struct
__device__ KQED_PRIVATE
struct invariants
set_invariants( const double xv[4] ,
		const double yv[4] ,
		const struct Grid_coeffs Grid ) ;

#endif
