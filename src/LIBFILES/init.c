/**
   @file init.c
   @brief a place where things are initialised
 */
#include "KQED.h"

#include "corr_malloc.h"
#include "io.h"        // read_ff()
#include "pi_pert.h"   // init_g8()

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

  return 0 ;
}
