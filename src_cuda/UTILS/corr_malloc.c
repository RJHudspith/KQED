/**
   @file corr_malloc.c
   @brief memory allocation wrapper
 */
#define _POSIX_C_SOURCE 200809L

#include "KQED.h"

// memalign wrapper
int 
corr_malloc( void **memptr , 
	     const size_t alignment , 
	     const size_t size )
{
#if (defined HAVE_IMMINTRIN_H)
  return posix_memalign( memptr , alignment , size ) ;
#else
  *memptr = malloc( size ) ;
  return ( *memptr == NULL ) ? 1 : 0 ;
#endif
}
