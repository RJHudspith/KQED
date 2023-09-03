#ifndef TABD_H
#define TABD_H

__device__
int
Tabd_xeq0( const double yv[4] ,
	   const struct Grid_coeffs Grid ,
	   double tI[4][4][4] ,
	   double tII[4][4][4] ,
	   double tIII[4][4][4] ) ;

__device__
int
Tabd_yeq0( const double xv[4] ,
	   const struct Grid_coeffs Grid ,
	   double tI[4][4][4] ,
	   double tII[4][4][4] ,
	   double tIII[4][4][4] ) ;

#endif
