/**
   @file chnr_dV.c
   @brief contains the function chnr_dV
 */
#include "KQED.h"

#include "getff-new.h" // accessv and extractff

// little struct definition to store all these variables
struct vtmp {
  double g2, dg1dx, dg1dy, dg1dcb, dg2dx, dg2dy, dg2dcb , dg3dx, dg3dy, dg3dcb;
  double ddg2dxdx, ddg2dxdy, ddg2dxdcb, ddg2dydcb, ddg2dcbdcb;
  double ddg3dxdx, ddg3dxdy, ddg3dxdcb, ddg3dydcb, ddg3dcbdcb;
  double ddg1dxdx, ddg1dxdy, ddg1dxdcb, ddg1dydcb, ddg1dcbdcb;
  double ddg1dydy , ddg2dydy ;
} ;

static int
taylorx_contrib( const struct Grid_coeffs Grid ,
		 struct vtmp *D ,
		 const double x ,
		 const double y ,
		 const double xsq ,
		 const bool flag2 ,
		 double cb )
{ 
  const int iy_tay  = bsrch( Grid.TX[ YY ] , y , 0 , Grid.NY_tay ) ;
  const int iy = bsrch( Grid.YY , y , 0 , Grid.nstpy ) ;
  
  if( iy_tay == Grid.NY_tay-1 || iy == Grid.nstpy-1 ) {
    fprintf( stderr , "chnr_dV y is above maximum grid value\n" ) ;
    fprintf( stderr , "%f >= %f \n" , y , Grid.TX[YY][iy_tay] ) ;
    fprintf( stderr , "%f >= %f \n" , y , Grid.YY[iy] ) ;
    return 1 ;
  }

  const int iy2_tay = iy_tay+1;
  const int ix=0;
  
  const double ay = (Grid.TX[ YY ][iy2_tay]-y)/(Grid.TX[YY][iy2_tay]-Grid.TX[YY][iy_tay]);
  const double ax = (Grid.XX[0]-x)/(Grid.XX[0]);

  const double xb=Grid.XX[0] ;
  const double xbsq=xb*xb;

  const double g2a = 0.0;
  const double g2b = accessv( false, true, ix, iy, QG2 , false, d0cb, cb, y, Grid );
  D -> g2 = lerp( ax , g2a , g2b ) ;

  const double ddg2adxdcb = lerp( ay , Grid.TX[ G21 ][iy_tay] , Grid.TX[ G21 ][iy2_tay] ) ;
  const double ddg2bdxdcb = accessv( false, true, ix, iy, dxQG2, false, d1cb, cb, y, Grid );
  D -> ddg2dxdcb = lerp( ax , ddg2adxdcb , ddg2bdxdcb ) ;

  const double dg2bdx = accessv( false, true, ix, iy, dxQG2, false, d0cb, cb, y , Grid );
  D -> dg2dx = D -> ddg2dxdcb*cb; // ax*dg2adx + (1.0-ax)*dg2bdx;

  const double dg2bdcb = accessv( false, true , ix, iy, QG2, false, d1cb, cb, y , Grid );
  D -> dg2dcb = D -> ddg2dxdcb*x; // dg2dx*x/cb; // ax*dg2adcb + (1.0-ax)*dg2bdcb;

  const double dg2bdy = accessv( false, false, ix, iy, QG2, true, d0cb, cb, y , Grid );
  D -> dg2dy = lerp( ax , 0.0 , dg2bdy ) ;

  const double ddg2adxdy = cb*lerp( ay , Grid.TX[ G21dy ][iy_tay] , Grid.TX[ G21dy ][iy2_tay] ) ;
  const double ddg2bdxdy = accessv( false, false, ix, iy, dxQG2, true, d0cb, cb, y , Grid );
  D -> ddg2dxdy = lerp( ax , ddg2adxdy , ddg2bdxdy ) ;

  const double ddg2adxdx = 2.0*( lerp( ay , Grid.TX[ G22A ][iy_tay] , Grid.TX[ G22A ][iy2_tay] ) 
				 + lerp( ay , Grid.TX[ G22B ][iy_tay] , Grid.TX[ G22B ][iy2_tay] )
				 *(0.4+0.6*(2.0*cb*cb-1)));
  const double ddg2bdxdx = accessv( false, true , ix, iy, d2xQG2, false, d0cb, cb, y , Grid );
  D -> ddg2dxdx = lerp( ax , ddg2adxdx , ddg2bdxdx ) ;


  const double ddg2bdcbdcb = accessv(false, true, ix, iy, QG2, false, d2cb, cb, y , Grid );
  D -> ddg2dcbdcb = ddg2bdcbdcb*xsq/xbsq; // the function starts quadratically

  const double ddg2bdydcb = accessv( false, false, ix, iy, QG2, true, d1cb, cb, y, Grid );
  D -> ddg2dydcb = lerp( ax , 0.0 , ddg2bdydcb ) ;
  
  const double dg3bdy      = accessv( false, false, ix, iy, QG3   , true , d0cb, cb, y , Grid );
  const double ddg3bdxdy   = accessv( false, false, ix, iy, dxQG3 , true , d0cb, cb, y , Grid );
  const double ddg3bdxdcb  = accessv( false, true , ix, iy, dxQG3 , false, d1cb, cb, y , Grid );
  const double ddg3bdxdx   = accessv( false, true , ix, iy, d2xQG3, false, d0cb, cb, y , Grid );
  const double ddg3bdcbdcb = accessv( false, true , ix, iy, QG3   , false, d2cb, cb, y , Grid );
  const double ddg3bdydcb  = accessv( false, false, ix, iy, QG3   , true , d1cb, cb, y , Grid );

  const double dg1bdy = dg3bdy - cb/xb*g2b -(y/xb)*cb*dg2bdy;
  const double ddg1bdxdx =  ddg3bdxdx -2.0*y/(xbsq*xb)*cb*g2b+y/xbsq*cb*dg2bdx + y/xbsq*cb*dg2bdx - y/xb*cb*ddg2bdxdx;
  const double ddg1bdxdy = ddg3bdxdy +cb/xbsq*g2b + y*cb/xbsq*dg2bdy -cb/xb*dg2bdx - y/xb*cb*ddg2bdxdy;
  const double ddg1bdxdcb = ddg3bdxdcb +y/xbsq*g2b +y/xbsq*cb*dg2bdcb -y/xb*dg2bdx -y/xb*cb*ddg2bdxdcb;
  const double ddg1bdydcb = ddg3bdydcb - g2b/xb -cb/xb*dg2bdcb -y/xb*dg2bdy - y/xb*cb*ddg2bdydcb;
  const double ddg1bdcbdcb = ddg3bdcbdcb - 2.0*y/xb*dg2bdcb - y/xb*cb*ddg2bdcbdcb;
    
  const double dg1ady = lerp( ay , Grid.TX[ G3Ady ][iy_tay] , Grid.TX[ G3Ady ][iy2_tay] )
    - lerp( ay , Grid.TX[ G3Bdy ][iy_tay] , Grid.TX[ G3Bdy ][iy2_tay] ) ;
  
  const double ddg1adxdy = cb*( lerp( ay , Grid.TX[ G31Ady ][iy_tay] , Grid.TX[ G31Ady ][iy2_tay] )
				-(2/3.)*lerp( ay , Grid.TX[ G31Bdy ][iy_tay] , Grid.TX[ G31Bdy ][iy2_tay] )
				-lerp( ay , Grid.TX[ G22A ][iy_tay] , Grid.TX[ G22A ][iy2_tay] )
				-y*lerp( ay , Grid.TX[ G22Ady ][iy_tay] , Grid.TX[ G22Ady ][iy2_tay] )) ;
  const double ddg1adxdcb =
    lerp( ay , Grid.TX[ G31A ][iy_tay] , Grid.TX[ G31A ][iy2_tay] )
    -(2/3.)*lerp( ay , Grid.TX[ G31B ][iy_tay] , Grid.TX[ G31B ][iy2_tay] )
    -y*lerp( ay , Grid.TX[ G22A ][iy_tay] , Grid.TX[ G22A ][iy2_tay] ) ;

  const double ddg1adxdx = ddg1bdxdx; 
  D -> ddg1dxdcb  = lerp( ax , ddg1adxdcb , ddg1bdxdcb ) ;
  D -> dg1dx      = D -> ddg1dxdcb*cb; // ax*dg1adx + (1.0-ax)*dg1bdx;     
  D -> dg1dcb     = D -> ddg1dxdcb*x; // ax*dg1adcb + (1.0-ax)*dg1bdcb;     
  D -> dg1dy      = lerp( ax , dg1ady , dg1bdy ) ;
  D -> ddg1dxdy   = lerp( ax , ddg1adxdy , ddg1bdxdy ) ;
  D -> ddg1dxdx   = lerp( ax , ddg1adxdx , ddg1bdxdx ) ;
  D -> ddg1dcbdcb = ddg1bdcbdcb*xsq/xbsq; // the function starts quadratically 
  D -> ddg1dydcb  = lerp( ax , 0.0 , ddg1bdydcb ) ;

#ifdef VERBOSE
  // the following don't get used but computed for regression to previous code
  const double dg2ady = 0.0; // (ay*G2dy_tay[iy_tay] + (1.0-ay)*G2dy_tay[iy2_tay])*cb;
  const double dg3bdx  = accessv( false, true , ix, iy, dxQG3 , false, d0cb, cb, y , Grid );
  const double dg3bdcb = accessv( false, true , ix, iy, QG3   , false, d1cb, cb, y , Grid );
  const double dg2adx = ddg2adxdcb*cb;
  const double dg1adx = cb*( lerp( ay , Grid.TX[ G31A ][iy_tay] , Grid.TX[ G31A][iy2_tay] )
			     -(2/3.)*lerp( ay , Grid.TX[ G31B ][iy_tay] , Grid.TX[ G31B ][iy2_tay] )
			     -y*lerp( ay , Grid.TX[ G22A ][iy_tay] , Grid.TX[ G22A ][iy2_tay] ) ) ;
  const double dg1bdcb = dg3bdcb - y/xb*g2b - y/xb*cb*dg2bdcb;
  const double dg1bdx = dg3bdx + y/xbsq*cb*g2b - y/xb*cb*dg2bdx;
  const double dg2adcb = 0.0; // (ay*G21_tay[iy_tay] + (1.0-ay)*G21_tay[iy2_tay])*x; 
  const double dg1adcb = 0.0;
  const double ddg1adcbdcb = 0.0;
  const double ddg2adcbdcb = 0.0;
  const double ddg2adydcb = 0.0;
  const double ddg1adydcb = 0.0;
  fprintf( stdout , "ddg2adxdcb= %.11lg  ddg2bdxdcb= %.11lg\n", ddg2adxdcb, ddg2bdxdcb);
  fprintf( stdout , "dg2adx= %.11lg  dg2bdx= %.11lg\n", dg2adx, dg2bdx); 
  fprintf( stdout , "dg2adcb= %.11lg  dg2bdcb= %.11lg\n", dg2adcb, dg2bdcb); 
  fprintf( stdout , "dg2ady= %.11lg  dg2bdy= %.11lg\n", dg2ady, dg2bdy); 
  fprintf( stdout , "ddg2adxdy= %.11lg  ddg2bdxdy= %.11lg\n", ddg2adxdy, ddg2bdxdy);
  fprintf( stdout , "ddg2adxdx= %.11lg  ddg2bdxdx= %.11lg\n", ddg2adxdx, ddg2bdxdx);
  fprintf( stdout , "ddg2adcbdcb= %.11lg  ddg2bdcbdcb= %.11lg\n", ddg2adcbdcb, ddg2bdcbdcb);
  fprintf( stdout , "ddg2adydcb= %.11lg  ddg2bdydcb= %.11lg\n", ddg2adydcb, ddg2bdydcb);
  fprintf( stdout , "dg1adx= %.11lg  dg1bdx= %.11lg\n", dg1adx, dg1bdx);
  fprintf( stdout , "dg1ady= %.11lg  dg1bdy= %.11lg\n", dg1ady, dg1bdy);
  fprintf( stdout , "dg1adcb= %.11lg  dg1bdcb= %.11lg\n", dg1adcb, dg1bdcb);
  fprintf( stdout , "ddg1adxdy= %.11lg ddg1bdxdy= %.11lg\n", ddg1adxdy, ddg1bdxdy);
  fprintf( stdout , "ddg1adxdcb= %.11lg ddg1bdxdcb= %.11lg\n", ddg1adxdcb, ddg1bdxdcb);
  fprintf( stdout , "ddg1adxdx= %.11lg ddg1bdxdx= %.11lg\n", ddg1adxdx, ddg1bdxdx);
  fprintf( stdout , "ddg1adcbdcb= %.11lg ddg1bdcbdcb= %.11lg\n", ddg1adcbdcb, ddg1bdcbdcb);
  fprintf( stdout , "ddg1adydcb= %.11lg ddg1bdydcb= %.11lg\n", ddg1adydcb, ddg1bdydcb);
#endif // VERBOSE

  if(flag2) {
    // in the case of a swap, you need the second derivative wrt y
    const double ya = Grid.YY[iy];
    const double yb = Grid.YY[iy+1] ; //ya+Grid.ystp;
    const double ddg2adydy= 0.0;

    const double dg2bdy_ya = accessv( false, false, ix, iy, QG2, true, d0cb, cb, ya, Grid );
    const double dg2bdy_yb = accessv( false, false, ix, iy, QG2, true, d0cb, cb, yb, Grid );
    const double ddg2bdydy = (dg2bdy_yb-dg2bdy_ya)/Grid.ystp;
    D -> ddg2dydy = ax*ddg2adydy + (1.0-ax)*ddg2bdydy;
       
    const double ddg1adydy = (Grid.TX[ G3Ady ][iy2_tay]- Grid.TX[ G3Bdy ][iy2_tay]-\
			      ( Grid.TX[ G3Ady ][iy_tay]- Grid.TX[ G3Bdy ][iy_tay]))/(Grid.TX[YY][iy2_tay] - Grid.TX[YY][iy_tay] ) ;

    const double dg3bdy_ya = accessv( false, false, ix, iy, QG3, true, d0cb, cb, ya, Grid );
    const double dg3bdy_yb = accessv( false, false, ix, iy, QG3, true, d0cb, cb, yb, Grid );

    const double ddg3bdydy = (dg3bdy_yb-dg3bdy_ya)/Grid.ystp;
    const double ddg1bdydy = ddg3bdydy -(y*ddg2bdydy+2.0*dg2bdy)*cb/x;
    D -> ddg1dydy = ax*ddg1adydy + (1.0-ax)*ddg1bdydy;
  }

  return 0 ; 
}

// flips the sign of all the dg2 values
static struct vtmp
v_flip2( const struct vtmp D )
{
  struct vtmp F = D ; // memcpy
  F.g2         = -D.g2 ;
  F.dg2dx      = -D.dg2dx ;
  F.dg2dy      = -D.dg2dy ;
  F.dg2dcb     = -D.dg2dcb ;
  F.ddg2dxdx   = -D.ddg2dxdx ;
  F.ddg2dxdy   = -D.ddg2dxdy ;
  F.ddg2dydy   = -D.ddg2dydy ;
  F.ddg2dxdcb  = -D.ddg2dxdcb ;
  F.ddg2dydcb  = -D.ddg2dydcb ;
  F.ddg2dcbdcb = -D.ddg2dcbdcb ;  
  return F ;
}

static void
XMYSWAP_V_contrib( struct vtmp *D ,
		   const double x ,
		   const double xmy ,
		   const double y ,
		   const double xsq ,
		   const double xmysq ,
		   const double ysq ,
		   const double ca ,
		   const double cborig )
{ 
  const struct vtmp F = v_flip2( *D ) ; 
  const double cd = (y-x*cborig)/xmy;
  
  D -> g2     = F.g2;
  D -> dg2dx  = F.dg2dx + (1.0-ca*ca)/xmy*F.dg2dcb + ca * F.dg2dy;
  D -> dg2dcb = -(ysq*cd*F.dg2dcb + x*y*xmy*F.dg2dy)/xmysq;
  D -> dg2dy  = -(cborig+ca*cd)/xmy*F.dg2dcb+cd*F.dg2dy;

  D -> dg1dx  = F.dg1dx + (1.0-ca*ca)/xmy*F.dg1dcb + ca * F.dg1dy - D -> dg2dx;
  D -> dg1dcb = -(ysq*cd*F.dg1dcb + x*y*xmy*F.dg1dy)/xmysq - D -> dg2dcb;
  D -> dg1dy  = -(cborig+ca*cd)/xmy*F.dg1dcb+cd*F.dg1dy - D -> dg2dy;

  // ddg2*
  D -> ddg2dxdx = (F.ddg2dxdx + (1.0-ca*ca)/xmy*F.ddg2dxdcb + ca * F.ddg2dxdy)
    - 3.0*ca*(1.0-ca*ca)/xmysq*F.dg2dcb
    + (1.0-ca*ca)/xmy*(F.ddg2dxdcb + (1.0-ca*ca)/xmy*F.ddg2dcbdcb + ca * F.ddg2dydcb)
    + (1.0-ca*ca)/xmy*F.dg2dy + ca*(F.ddg2dxdy + (1.0-ca*ca)/xmy*F.ddg2dydcb + ca * F.ddg2dydy);
     
  D -> ddg2dxdcb = ysq/(xmysq*xmy)*(cborig+3.0*ca*cd)*F.dg2dcb
    -ysq/xmysq*cd*(F.ddg2dxdcb + (1.0-ca*ca)/xmy*F.ddg2dcbdcb + ca * F.ddg2dydcb)
    + y/xmysq*(x*ca-xmy)*F.dg2dy-x*y/xmy*(F.ddg2dxdy + (1.0-ca*ca)/xmy*F.ddg2dydcb
					  + ca * F.ddg2dydy);
  
  D -> ddg2dxdy = (-(cborig+ca*cd)/xmy*F.ddg2dxdcb+cd*F.ddg2dxdy)
    +(2.0*ca*cborig+cd*(3.0*ca*ca-1.0))/xmysq*F.dg2dcb
    +(1.0 - ca*ca)/xmy*(-(cborig + ca*cd)/xmy*F.ddg2dcbdcb + cd*F.ddg2dydcb)
    -(cborig+ca*cd)/xmy*F.dg2dy + ca*(-(cborig+ca*cd)/xmy*F.ddg2dydcb+cd*F.ddg2dydy);

  D -> ddg2dydcb = y/(xmysq*xmy)*((3.0*cd*cd-1.0)*y-2*xmy*cd)*F.dg2dcb
    -ysq/(xmysq)*cd*(-(cborig+ca*cd)/xmy*F.ddg2dcbdcb+cd*F.ddg2dydcb)
    + x/xmysq*(y*cd-xmy)*F.dg2dy
    -x*y/xmy*(-(cborig+ca*cd)/xmy*F.ddg2dydcb+cd*F.ddg2dydy);

  D -> ddg2dcbdcb = x*ysq/(xmysq*xmysq)*(xmy-3.0*y*cd)*F.dg2dcb
    - ysq/(xmysq)*cd*(-(ysq*cd*F.ddg2dcbdcb + x*y*xmy*F.ddg2dydcb)/xmysq)
    - xsq*ysq/(xmysq*xmy)*F.dg2dy
    - x*y/xmy*(-(ysq*cd*F.ddg2dydcb + x*y*xmy*F.ddg2dydy)/xmysq);

  // these have the 1s in
  D -> ddg1dxdx = ( F.ddg1dxdx + (1.0-ca*ca)/xmy*F.ddg1dxdcb + ca * F.ddg1dxdy)
    - 3.0*ca*(1.0-ca*ca)/xmysq*F.dg1dcb
    + (1.0-ca*ca)/xmy*(F.ddg1dxdcb + (1.0-ca*ca)/xmy*F.ddg1dcbdcb + ca * F.ddg1dydcb)
    + (1.0-ca*ca)/xmy*F.dg1dy + ca*(F.ddg1dxdy + (1.0-ca*ca)/xmy*F.ddg1dydcb + ca * F.ddg1dydy)
    - D -> ddg2dxdx;
  
  D -> ddg1dxdcb = ysq/(xmysq*xmy)*(cborig+3.0*ca*cd)*F.dg1dcb
    - ysq/xmysq*cd*(F.ddg1dxdcb + (1.0-ca*ca)/xmy*F.ddg1dcbdcb + ca * F.ddg1dydcb)
    + y/xmysq*(x*ca-xmy)*F.dg1dy
    - x*y/xmy*(F.ddg1dxdy+ (1.0-ca*ca)/xmy*F.ddg1dydcb + ca * F.ddg1dydy)
    - D -> ddg2dxdcb;

  D -> ddg1dxdy = (-(cborig+ca*cd)/xmy*F.ddg1dxdcb+cd*F.ddg1dxdy)
    +(2.0*ca*cborig+cd*(3.0*ca*ca-1.0))/xmysq*F.dg1dcb
    +(1.0 - ca*ca)/xmy*(-(cborig + ca*cd)/xmy*F.ddg1dcbdcb + cd*F.ddg1dydcb)
    -(cborig+ca*cd)/xmy*F.dg1dy + ca*(-(cborig+ca*cd)/xmy*F.ddg1dydcb+cd*F.ddg1dydy)
    - D -> ddg2dxdy;

  D -> ddg1dydcb = y/(xmysq*xmy)*((3.0*cd*cd-1.0)*y-2*xmy*cd)*F.dg1dcb
    -ysq/(xmysq)*cd*(-(cborig+ca*cd)/xmy*F.ddg1dcbdcb+cd*F.ddg1dydcb)
    + x/xmysq*(y*cd-xmy)*F.dg1dy-x*y/xmy*(-(cborig+ca*cd)/xmy*F.ddg1dydcb+cd*F.ddg1dydy)
    - D -> ddg2dydcb;

  D -> ddg1dcbdcb = x*ysq/(xmysq*xmysq)*(xmy-3.0*y*cd)*F.dg1dcb
    - ysq/(xmysq)*cd*(-(ysq*cd*F.ddg1dcbdcb + x*y*xmy*F.ddg1dydcb)/xmysq)
    - xsq*ysq/(xmysq*xmy)*F.dg1dy
    - x*y/xmy*(-(ysq*cd*F.ddg1dydcb + x*y*xmy*F.ddg1dydy)/xmysq)
    - D -> ddg2dcbdcb;
  
  return ;
}

// using the computed values in D sets the array dv
static void
get_dv( const double x ,
	const double xv[4] ,
	const double y ,
	const double yv[4] ,
	const double cb ,
	struct vtmp D ,
	double dv[4][4][4] )
{
  // look up tables and general precomputations
  const double xhat[4] = { xv[0]/x , xv[1]/x , xv[2]/x , xv[3]/x } ;
  const double yhat[4] = { yv[0]/y , yv[1]/y , yv[2]/y , yv[3]/y } ;
  const double c2v[4] = { (xhat[0]-cb*yhat[0])/y , (xhat[1]-cb*yhat[1])/y ,
			  (xhat[2]-cb*yhat[2])/y , (xhat[3]-cb*yhat[3])/y } ;
  const double c4v[4] = { (yhat[0]-cb*xhat[0])/x , (yhat[1]-cb*xhat[1])/x ,
			  (yhat[2]-cb*xhat[2])/x , (yhat[3]-cb*xhat[3])/x } ;
  const double D1[4] = { xv[0]*D.ddg1dxdx+yv[0]*D.ddg2dxdx ,
			 xv[1]*D.ddg1dxdx+yv[1]*D.ddg2dxdx ,
			 xv[2]*D.ddg1dxdx+yv[2]*D.ddg2dxdx ,
			 xv[3]*D.ddg1dxdx+yv[3]*D.ddg2dxdx } ;
  const double D2[4] = { xv[0]*D.ddg1dxdy+yv[0]*D.ddg2dxdy ,
			 xv[1]*D.ddg1dxdy+yv[1]*D.ddg2dxdy ,
			 xv[2]*D.ddg1dxdy+yv[2]*D.ddg2dxdy ,
			 xv[3]*D.ddg1dxdy+yv[3]*D.ddg2dxdy } ;
  const double D3[4] = { xv[0]*D.ddg1dxdcb+yv[0]*D.ddg2dxdcb ,
			 xv[1]*D.ddg1dxdcb+yv[1]*D.ddg2dxdcb ,
			 xv[2]*D.ddg1dxdcb+yv[2]*D.ddg2dxdcb ,
			 xv[3]*D.ddg1dxdcb+yv[3]*D.ddg2dxdcb } ;
  const double D4[4] = { xv[0]*D.ddg1dydcb+yv[0]*D.ddg2dydcb ,
			 xv[1]*D.ddg1dydcb+yv[1]*D.ddg2dydcb ,
			 xv[2]*D.ddg1dydcb+yv[2]*D.ddg2dydcb ,
			 xv[3]*D.ddg1dydcb+yv[3]*D.ddg2dydcb } ;
  const double D5[4] = { xv[0]*D.ddg1dxdcb+yv[0]*D.ddg2dxdcb ,
			 xv[1]*D.ddg1dxdcb+yv[1]*D.ddg2dxdcb ,
			 xv[2]*D.ddg1dxdcb+yv[2]*D.ddg2dxdcb ,
			 xv[3]*D.ddg1dxdcb+yv[3]*D.ddg2dxdcb } ;
  const double D6[4] = { xv[0]*D.ddg1dcbdcb+yv[0]*D.ddg2dcbdcb ,
			 xv[1]*D.ddg1dcbdcb+yv[1]*D.ddg2dcbdcb ,
			 xv[2]*D.ddg1dcbdcb+yv[2]*D.ddg2dcbdcb ,
			 xv[3]*D.ddg1dcbdcb+yv[3]*D.ddg2dcbdcb } ;
  const double D7[4] = { (xv[0]*D.dg1dx+yv[0]*D.dg2dx)/x ,
			 (xv[1]*D.dg1dx+yv[1]*D.dg2dx)/x ,
			 (xv[2]*D.dg1dx+yv[2]*D.dg2dx)/x ,
			 (xv[3]*D.dg1dx+yv[3]*D.dg2dx)/x } ;
  const double D8[4] = { (xv[0]*D.dg1dcb+yv[0]*D.dg2dcb)/x ,
			 (xv[1]*D.dg1dcb+yv[1]*D.dg2dcb)/x ,
			 (xv[2]*D.dg1dcb+yv[2]*D.dg2dcb)/x ,
			 (xv[3]*D.dg1dcb+yv[3]*D.dg2dcb)/x } ;

  int bet , alf , dta ;
  for(alf=0;alf<4;alf++)  {
    const double astuff = (D.dg1dx+D.dg2dx)*xhat[alf]+(D.dg1dcb+D.dg2dcb)*c4v[alf] ;
    for(bet=0;bet<4;bet++)  {
      // this only depends on beta
      const double bstuff = (D.dg1dx*xhat[bet]+D.dg1dy*yhat[bet]
			     +D.dg1dcb*(xhat[bet]/y+yhat[bet]/x
					-cb*(xhat[bet]/x+yhat[bet]/y))) ;
      const double c2c4 = c2v[bet]+c4v[bet] ;
      const int d_ab = ( alf==bet ) ;
      const double xaxb = xhat[alf]*xhat[bet] ;
      const double yayb = yhat[alf]*yhat[bet] ;
      const double xayb = xhat[alf]*yhat[bet] ;
      const double yaxb = yhat[alf]*xhat[bet] ;
      const double c1 = ((cb*(3.0*xaxb-d_ab)-xayb-yaxb)/x+(d_ab+xayb*cb-xaxb-yayb)/y) ;
      
      for(dta=0;dta<4;dta++) {
	const int d_bd = (bet==dta) ;
	const int d_ad = (alf==dta) ;
	
	register double f = d_bd*astuff; 
	 
	f += d_ad*bstuff;
	f += D1[dta]*xaxb+D2[dta]*xayb;
	f += (D3[dta]*xhat[bet]+D4[dta]*yhat[bet])*c4v[alf];
	f += (D5[dta]*xhat[alf]+D6[dta]*c4v[alf])*c2c4 ;
	f += D7[dta]*(d_ab-xaxb);

	dv[alf][bet][dta] = f + D8[dta]*c1 ;
      }
    }
  }
  return ;
}

// returns dv[alf][bet][dta] = \partial^(x)_\alpha (\partial^(x)_\beta + \partial^(y)_\beta) < \epsilon_\delta I>_\epsilon
int
chnr_dV( const double xv[4] ,
	 const double yv[4] ,
	 const struct invariants Inv ,
	 const struct Grid_coeffs Grid ,
	 const struct AVX_precomps *PC ,
	 double dv[4][4][4] )
{
  struct vtmp D = {} ; // zero this guy

  if( Inv.x < Grid.XX[0] ) {
    // if this messes up we leave
    if( taylorx_contrib( Grid , &D , Inv.x , Inv.y , Inv.xsq ,
			 Inv.flag2 , Inv.cb ) == 1 ) {
      return 1 ;
    }
  } else if( Inv.x > Grid.XX[ Grid.nstpx -1 ] ) {
    return 1 ;
  } else {

    // this one is all on its own and that is sad
    D.dg3dy = extractff( QG3   , true , d0cb, Inv, Grid , PC );
    
    // use the new extract code
    double f[4] KQED_ALIGN ;
    extractff2( QG2 , d0cb , Inv , Grid , PC , f ) ;
    D.g2    = f[3] ;
    D.dg2dy = f[2] ;

    extractff2( dxQG2 , d0cb , Inv , Grid , PC , f ) ;
    D.ddg2dxdx  = f[1] ;
    D.dg2dx     = f[3] ;
    D.ddg2dxdy  = f[2] ;
    D.ddg2dxdcb = f[0] ;

    extractff2( dxQG3 , d0cb , Inv , Grid , PC , f ) ;
    D.ddg3dxdx  = f[1] ;
    D.dg3dx     = f[3] ;
    D.ddg3dxdy  = f[2] ;
    D.ddg3dxdcb = f[0] ;
    
    extractff2( QG2 , d1cb , Inv , Grid , PC , f ) ;
    D.dg2dcb     = f[3] ;
    D.ddg2dydcb  = f[2] ;
    D.ddg2dcbdcb = f[0] ;

    extractff2( QG3 , d1cb , Inv , Grid , PC , f ) ;
    D.dg3dcb     = f[3] ;
    D.ddg3dydcb  = f[2] ;
    D.ddg3dcbdcb = f[0] ;

    D.dg1dx  = D.dg3dx + Inv.y/Inv.xsq*Inv.cb*D.g2 - Inv.y/Inv.x*Inv.cb*D.dg2dx;
    D.dg1dy  = D.dg3dy - Inv.cb/Inv.x*D.g2 -(Inv.y/Inv.x)*Inv.cb*D.dg2dy;
    D.dg1dcb = D.dg3dcb - Inv.y/Inv.x*D.g2 - Inv.y/Inv.x*Inv.cb*D.dg2dcb;
    
    D.ddg1dxdx   =  D.ddg3dxdx
      - 2.0*Inv.y/(Inv.xsq*Inv.x)*Inv.cb*D.g2
      + Inv.y/Inv.xsq*Inv.cb*D.dg2dx
      + Inv.y/Inv.xsq*Inv.cb*D.dg2dx
      - Inv.y/Inv.x*Inv.cb*D.ddg2dxdx;

    D.ddg1dxdy   = D.ddg3dxdy +Inv.cb/Inv.xsq*D.g2
      + Inv.y*Inv.cb/Inv.xsq*D.dg2dy
      - Inv.cb/Inv.x*D.dg2dx
      - Inv.y/Inv.x*Inv.cb*D.ddg2dxdy;

    D.ddg1dxdcb  = D.ddg3dxdcb
      + Inv.y/Inv.xsq*D.g2
      + Inv.y/Inv.xsq*Inv.cb*D.dg2dcb
      - Inv.y/Inv.x*D.dg2dx
      - Inv.y/Inv.x*Inv.cb*D.ddg2dxdcb;
    
    D.ddg1dydcb  = D.ddg3dydcb
      - D.g2/Inv.x
      - Inv.cb/Inv.x*D.dg2dcb
      - Inv.y/Inv.x*D.dg2dy
      - Inv.y/Inv.x*Inv.cb*D.ddg2dydcb;
    
    D.ddg1dcbdcb = D.ddg3dcbdcb
      - 2.0*Inv.y/Inv.x*D.dg2dcb
      - Inv.y/Inv.x*Inv.cb*D.ddg2dcbdcb;
      
    if( Inv.flag2 ) {
      // in case of swap, you need the second derivative wrt y
      const int iy1 = bsrch( Grid.YY , Inv.y , 0 , Grid.nstpy ) ;

      if( iy1 == Grid.nstpy-1 ) {
	fprintf( stderr , "chnr_dV xmyswap y is at or above edge of grid\n" ) ;
	fprintf( stderr , "%f > %f \n" , Inv.y , Grid.YY[ Grid.nstpy-1] ) ;
	return 1 ;
      }
      
      struct invariants Inv1 = Inv , Inv2 = Inv ;
      Inv1.y = Grid.YY[iy1] ;
      Inv2.y = Grid.YY[iy1+1] ;
      
      double fq1 = extractff( QG2, true, d0cb, Inv1, Grid, PC );
      double fq2 = extractff( QG2, true, d0cb, Inv2, Grid, PC );
      D.ddg2dydy = (fq2-fq1)/Grid.ystp;
      
      fq1 = extractff( QG3, true, d0cb, Inv1, Grid , PC );
      fq2 = extractff( QG3, true, d0cb, Inv2, Grid , PC );
      const double ddg3dydy = (fq2-fq1)/Grid.ystp; 
      D.ddg1dydy = ddg3dydy -(Inv.y*D.ddg2dydy+2.0*D.dg2dy)*Inv.cb/Inv.x;
    }
  }
  
  if( Inv.flag2) {
    // here, convert the derivatives of g2(x,ca,xmy) into the derivatives of g2(x,cb,y)
    XMYSWAP_V_contrib( &D , Inv.x , Inv.xmy , Inv.yorig ,
		       Inv.xsq , Inv.xmysq , Inv.ysq ,
		       Inv.cb , Inv.cborig ) ;
  } 
  
  get_dv( Inv.x , xv , Inv.yorig , yv , Inv.cborig , D , dv ) ;

  return 0 ;
}

