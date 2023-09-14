/**
   @file simd.h
   @brief various AVX/FMA SIMD function declarations for optimising the kernel
 */
#ifndef SIMD_H
#define SIMD_H

void
precompute_INV( struct intprecomp *INVy,
		const double y ,
		const double y1 ,
		const double y2 ,
		const size_t idx ) ;

#if (defined HAVE_IMMINTRIN_H) && (defined __AVX__)

void
accessv6( const int ix,
	  const FFidx nm,
	  const bool ndy,
	  const NDCB ndcb,
	  const double cb,
	  const struct Grid_coeffs Grid,
	  const struct AVX_precomps *PC,
	  const struct intprecomp INVy,
	  double f[4] ) ;

void
accessv7( const int ix,
	  const FFidx nm,
	  const NDCB ndcb,
	  const double cb,
	  const struct Grid_coeffs Grid,
	  const struct AVX_precomps *PC,
	  const struct intprecomp INVy,
	  double f[12] ) ;

int
init_LUT12( const size_t max_LUT ) ;

int
init_LUT34( const size_t max_LUT ) ;

int
init_LUT56( const size_t max_LUT ) ;

void
free_LUTs( void ) ;

#endif

#endif
