#ifndef CHNR_DT_H
#define CHNR_DT_H

// set dyv[alf][bet][dta]= \partial^(y)_\beta <(epsilon_alf epsilon_dta - 1/4 delta_{alf dta}) I>_epsilon
// and dxv[alf][bet][dta]= \partial^(x)_\beta <(epsilon_alf epsilon_dta - 1/4 delta_{alf dta}) I>_epsilon
int
chnr_dT( const double xv[4] ,
	 const double yv[4] ,
	 const struct invariants Inv ,
	 const struct Grid_coeffs Grid ,
	 const struct AVX_precomps *PC ,
	 double dxv[4][4][4] ,
	 double dyv[4][4][4] ) ;

#endif
