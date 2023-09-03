#ifndef CHNR_DS_H
#define CHNR_DS_H

// sets dyv[bet]= \partial^(y)_\beta < I >_epsilon
// sets dxv[bet]= \partial^(x)_\beta < I >_epsilon
__device__
int
chnr_dS( const double xv[4] ,
	 const double yv[4] ,
	 const struct invariants Inv ,
	 const struct Grid_coeffs Grid ,
	 double dxv[4] ,
	 double dyv[4] ) ;

#endif
