/**
   @file init.c
   @brief a place where things are initialised
 */
#include "KQED.h"

#include "corr_malloc.h"
#include "io.h"        // read_ff()
#include "pi_pert.h"   // init_g8()

#if (defined HAVE_IMMINTRIN_H) && (defined __AVX__)
  #include "simd.h" // allocation of LUTS
#endif

// set the Grid struct and all the grid information e.g. XX,YY
static int
set_grid( struct Grid_coeffs *Grid )
{  
  if( read_ff( Grid ) ) {
    return 1 ;
  }
  if( read_TAYLORX( Grid ) ) {
    return 1 ;
  }
  if( read_TAYLORY( Grid ) ) {
    return 1 ;
  }

#ifdef VERBOSE
  fprintf( stdout , "Successfully Initialized QED kernel FFs...\n");
#endif
  
  return 0 ;
}

// free the struct t
void
free_QED_temps( struct QED_kernel_temps *t )
{
  size_t a , i , j ;

#if (defined HAVE_IMMINTRIN_H) && (defined __AVX__)
  free_LUTs() ;
#endif

  // free these grid parameters
  if( t->Grid.XX != NULL ) {
    free( t->Grid.XX ) ;
  }
  if( t->Grid.nfx != NULL ) {
    free( t->Grid.nfx ) ;
  }
  if( t[0].Grid.YY != NULL ) {
    free( t->Grid.YY ) ;
  }

  // Free Ffm
  if( t->Grid.Ffm != NULL ) {
    for( a = 0 ; a < (size_t)t->Grid.Nffa ; a++ ) {
      if( t->Grid.Ffm[a] != NULL ) {
	for( i = 0 ; i < (size_t)t->Grid.nstpx ; i++ ) {
	  if( t->Grid.Ffm[a][i] != NULL ) {
	    for( j = 0 ; j < (size_t)t->Grid.nstpy ; j++ ) {
	      if( t->Grid.Ffm[a][i][j] != NULL ) {
		free( t->Grid.Ffm[a][i][j] ) ;
	      }
	    }
	    free( t->Grid.Ffm[a][i] ) ;
	  }
	}
	free( t->Grid.Ffm[a] ) ;
      }
    }
    free( t->Grid.Ffm ) ;
  }

  // Free Ffp
  if( t->Grid.Ffp != NULL ) {
    for( a = 0 ; a < (size_t)t->Grid.Nffa ; a++ ) {
      if( t->Grid.Ffp[a] != NULL ) {
	for( i = 0 ; i < (size_t)t->Grid.nstpx ; i++ ) {
	  if( t->Grid.Ffp[a][i] != NULL ) {
	    for( j = 0 ; j < (size_t)t->Grid.nstpy ; j++ ) {
	      if( t->Grid.Ffp[a][i][j] != NULL ) {
		free( t->Grid.Ffp[a][i][j] ) ;
	      }
	    }
	    free( t->Grid.Ffp[a][i] ) ;
	  }
	}
	free( t->Grid.Ffp[a] ) ;
      }
    }
    free( t->Grid.Ffp ) ;
  }
  
  // free the Taylor coefficients
  if( t->Grid.TX != NULL ) {
    for( i = 0 ; i < (size_t)t->Grid.NtayY ; i++ ) {
      if( t->Grid.TX[i] != NULL ) {
	free( t->Grid.TX[i] ) ;
      }
    }
    free( t->Grid.TX ) ;
  }
  
  if( t->Grid.TY != NULL ) {
    for( i = 0 ; i < (size_t)t->Grid.NtayX ; i++ ) {
      if( t->Grid.TY[i] != NULL ) {
	free( t->Grid.TY[i] ) ;
      }
    }
    free( t->Grid.TY ) ;
  }

  // clean up the allocations
  #if (defined HAVE_IMMINTRIN_H) && (defined __AVX__)
  if( t->Grid.PC != NULL ) {
    const size_t Nthreads = omp_get_max_threads() ;
    for( i = 0 ; i < t->Grid.Nffa*Nthreads ; i++ ) {
      if( t->Grid.PC[i].Fm != NULL ) {
	free( t->Grid.PC[i].Fm ) ;
      }
      if( t->Grid.PC[i].Fp != NULL ) {
	free( t->Grid.PC[i].Fp ) ;
      }     
    }
    free( t->Grid.PC ) ;
  }
  #endif
  
  // free G8
  if( t->G8 != NULL ) {
    free( t->G8 ) ;
  }
  
  return ;
}

// returns 1 if it messes up, returns 0 otherwise
int
initialise( struct QED_kernel_temps *t )
{
  // set the grid struct and initialise the ff arrays
  if( set_grid( &( t -> Grid ) ) ) {
    return 1 ;
  }

  t -> G8 = malloc( 65536 * sizeof( int ) ) ;
  init_g8( t );

#ifdef VERBOSE
  fprintf( stdout , "|INIT| Set g-factors\n") ;
#endif

  // initialise LUTs
#if (defined HAVE_IMMINTRIN_H) && (defined __AVX__)
  // assumes nfx is ordered, which it is otherwise we would
  // need to generically call a max
  const size_t max_LUT = t -> Grid.nfx[ t->Grid.nstpx-1 ] + 1 ;
  if( init_LUT12( max_LUT ) == 1 ) return 1 ;
  if( init_LUT34( max_LUT ) == 1 ) return 1 ;
  if( init_LUT56( max_LUT ) == 1 ) return 1 ;
  #ifdef VERBOSE
  fprintf( stdout , "LUTS initialised\n" ) ;
  #endif

  // init PC
  // find largest nfx
  size_t i , nfmax = 0 ;
  for( i = 0 ; i < (size_t)t->Grid.nstpx ; i++ ) {
    if( (size_t)t->Grid.nfx[i] > nfmax ) {
      nfmax = t->Grid.nfx[i] ;
    }
  }

  // allocate PC here instead of in the kernel computation
  const size_t Nthreads = (size_t)omp_get_max_threads() ;
  t->Grid.PC = malloc( t->Grid.Nffa*Nthreads*sizeof(struct AVX_precomps) ) ;  
  for( i = 0 ; i < t->Grid.Nffa*Nthreads ; i++ ) {
    if( corr_malloc( (void**)&t->Grid.PC[i].Fm , 32 , nfmax*sizeof(__m256d) ) != 0 ||
	corr_malloc( (void**)&t->Grid.PC[i].Fp , 32 , nfmax*sizeof(__m256d) ) != 0 ) {
      fprintf( stderr , "[INIT] PC allocation failed\n" ) ;
      return 1 ;
    }
  }
#endif
  
  return 0 ;
}
