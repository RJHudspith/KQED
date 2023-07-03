#ifndef CHNR_DV_H
#define CHNR_DV_H

// dv[alf][bet][dta] = \partial^(x)_\alpha (\partial^(x)_\beta + \partial^(y)_\beta) < \epsilon_\delta I>_\epsilon
__device__
int
chnr_dV( const double xv[4] ,
	 const double yv[4] ,
	 const struct invariants Inv ,
	 const struct Grid_coeffs Grid ,
	 double dv[4][4][4] ) ;

#endif
