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

#endif
