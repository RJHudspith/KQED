/**
   @file init.c
   @brief a place where things are initialised
 */
#include "KQED.h"

#include "io.h"        // read_ff()
#include "pi_pert.h"   // init_g8()

// set the Grid struct and all the grid information e.g. XX,YY.
// this includes uploading grid info to GPU, though Grid is still a host struct
// wrapping the device ptrs.
__host__
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
  // free these grid parameters
  if( t->Grid.XX != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.XX )) ;
  }
  if( t->Grid.nfx != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.nfx )) ;
  }
  if( t[0].Grid.YY != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.YY )) ;
  }

  // Free Ffm
  if( t->Grid.Ffm != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.Ffm )) ;
  }

  // Free Ffp
  if( t->Grid.Ffp != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.Ffp )) ;
  }
  
  // free the Taylor coefficients
  if( t->Grid.TX != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.TX )) ;
  }
  
  if( t->Grid.TY != NULL ) {
    checkCudaErrors(cudaFree( t->Grid.TY )) ;
  }

  // free G8
  if( t->G8 != NULL ) {
    checkCudaErrors(cudaFree( t->G8 )) ;
  }
  
  return ;
}

// returns 1 if it messes up, returns 0 otherwise
__host__
int
initialise( struct QED_kernel_temps *t )
{
  // set the grid struct and initialise the ff arrays
  if( set_grid( &( t -> Grid ) ) ) {
    return 1 ;
  }

  checkCudaErrors(cudaMalloc((void**)&t->G8, 65536 * sizeof( int ) )) ;
  init_g8( t );

#ifdef VERBOSE
  fprintf( stdout , "|INIT| Set g-factors\n") ;
#endif

  return 0 ;
}
