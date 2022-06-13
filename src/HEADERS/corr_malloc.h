/**
   @file corr_malloc.h
   @brief prototype declarations for malloc wrappers
 */
#ifndef CORR_MALLOC_H
#define CORR_MALLOC_H

/**
   @fn int corr_malloc( void **memptr , const size_t alignment , const size_t size )
   @brief aligned malloc wrapper
   defaults to malloc if HAVE_IMMINTRIN_H is not set
 */
int 
corr_malloc( void **memptr , 
	     const size_t alignment , 
	     const size_t size ) ;

#endif
