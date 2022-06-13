/**
   @file chnr_dS.c
   @brief contains the function chnr_dS
*/
#include "KQED.h"

#include "getff-new.h" // accessv and extractff

static int
compute_dg_TAYLORX( double *dg0dy ,
		    double *dg0dx ,
		    double *dg0dcb ,
		    const struct invariants Inv ,
		    const struct Grid_coeffs Grid )
{  
  const int iy_tay = bsrch( Grid.TX[ YY ] , Inv.y , 0 , Grid.NY_tay ) ;
  const int iy = bsrch( Grid.YY , Inv.y , 0 , Grid.nstpy ) ;
  
  if( iy_tay == Grid.NY_tay-1 || iy == Grid.nstpy-1 ) {
    fprintf( stderr , "chnr_dS y is at the end of the grid\n" ) ;
    fprintf( stderr , "%f >= %f \n" , Inv.y , Grid.TX[YY][iy_tay] ) ;
    fprintf( stderr , "%f >= %f \n" , Inv.y , Grid.YY[iy] ) ;
    return 1 ;
  }

  const int iy2_tay = iy_tay+1;
  const int ix  = 0;
  
  const double ay = (Grid.TX[YY][iy2_tay]-Inv.y)
    /(Grid.TX[YY][iy2_tay]-Grid.TX[YY][iy_tay]);
  const double ax = ( Grid.XX[0]-Inv.x)/( Grid.XX[0]);
  
  double f1 = Inv.cb*lerp( ay , Grid.TX[G0dx][iy_tay] , Grid.TX[G0dx][iy2_tay] ) ;  
  double f2 = accessv( false, true, ix, iy, dxQG0, false, d0cb, Inv.cb, Inv.y, Grid );
  *dg0dx = lerp( ax , f1 , f2 ) ;

  f1 = lerp( ay , Grid.TX[ G0dy ][iy_tay] , Grid.TX[ G0dy ][iy2_tay ] ) ;
  f2 = accessv( false, false, ix, iy, QG0, true , d0cb, Inv.cb, Inv.y, Grid );
  *dg0dy = lerp( ax , f1 , f2 ) ;
    
  f2 = accessv( false, true , ix, iy, QG0, false, d1cb, Inv.cb, Inv.y, Grid );
  *dg0dcb = lerp( ax , 0.0 , f2 ) ;

  return 0 ;
}

// sets dyv[bet]= \partial^(y)_\beta < I >_epsilon
// sets dxv[bet]= \partial^(x)_\beta < I >_epsilon
int
chnr_dS( const double xv[4] ,
	 const double yv[4] ,
	 const struct invariants Inv ,
	 const struct Grid_coeffs Grid ,
	 const struct AVX_precomps *PC ,
	 double dxv[4] ,
	 double dyv[4] )
{
  double dg0dy, dg0dx, dg0dcb;
  
  if( Inv.x < Grid.XX[0] ) {
    if( compute_dg_TAYLORX( &dg0dy , &dg0dx , &dg0dcb , Inv, Grid ) == 1 ) {
      return 1 ;
    }
  } else if( Inv.x > Grid.XX[ Grid.nstpx -1 ] ) {
    return 1 ;
  } else { 
    double f[4] KQED_ALIGN ;
    extractff2( QG0 , d0cb, Inv, Grid , PC , f );
    dg0dcb = f[0] ;
    dg0dx  = f[1] ;
    dg0dy  = f[2] ;
  }
  
  if( Inv.flag2 ) {
    const double dgs0dx  = dg0dx;
    const double dgs0dcb = dg0dcb;
    const double dgs0dy  = dg0dy;
    const double ca = Inv.cb; 
    const double cd = (Inv.yorig-Inv.x*Inv.cborig)/Inv.xmy;
    dg0dx  = dgs0dx + (1.0-ca*ca)/Inv.xmy*dgs0dcb + ca*dgs0dy;
    dg0dcb = -(Inv.ysq*cd*dgs0dcb + Inv.x*Inv.yorig*Inv.xmy*dgs0dy)/Inv.xmysq; 
    dg0dy  = -(Inv.cborig+ca*cd)/Inv.xmy*dgs0dcb+cd*dgs0dy;
  } 

  int bet ;
  for(bet=0;bet<4;bet++) {
    const double yhat = yv[bet]/Inv.yorig;
    const double xhat = xv[bet]/Inv.x;
    const double c2 = (xhat-Inv.cborig*yhat)/Inv.yorig;
    const double c4 = (yhat-Inv.cborig*xhat)/Inv.x;
    dyv[bet] = yhat*dg0dy+c2*dg0dcb;
    dxv[bet] = xhat*dg0dx+c4*dg0dcb;
  }
  return 0 ;
}
