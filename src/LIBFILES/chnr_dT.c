/**
   @file chnr_dT.c
   @brief contains the function chnr_dT
*/
#include "KQED.h"       // various definitions

#include "getff-new.h"  // accessv and extractff

// little storage for the temporary variables
struct ttmps {
  double ell1, ell2, ell3, v1, ell4, dv1dx, dv1dy, dv1dcb, dell2dx, dell2dy, dell2dcb;
  double dell1dx, dell1dy, dell1dcb, dell3dx, dell3dy, dell3dcb, dell4dy, dell4dx, dell4dcb;
} ;

static int
TAYLORX_contribution( struct ttmps *T ,
		      const struct invariants Inv ,
		      const struct Grid_coeffs Grid )
{
  // Just like non taylor-expanded version need to set this due to
  // it being set differently compared to the V and S
  const double ysq = Inv.flag2 ? Inv.xmysq : Inv.ysq ;
  
  // set the indices
  const int iy_tay  = bsrch( Grid.TX[ YY ] , Inv.y , 0 , Grid.NY_tay ) ;      
  const int iy = bsrch( Grid.YY , Inv.y , 0 , Grid.nstpy ) ;

  // check they are reasonable
  if( iy_tay == Grid.NY_tay-1 || iy == Grid.nstpy-1 ) {
    fprintf( stderr , "chnr_dT TAYLORX y is at or above edge of grid\n" ) ;
    fprintf( stderr , "%f >= %f\n" , Inv.y , Grid.TX[YY][iy_tay] ) ;
    fprintf( stderr , "%f >= %f\n" , Inv.y , Grid.YY[Grid.nstpy-1] ) ;
    return 1 ;
  }

  const int iy2_tay = iy_tay+1;
  const int ix = 0 ;
  
  const double ay = (Grid.TX[ YY ][iy2_tay]-Inv.y)/
    (Grid.TX[ YY ][iy2_tay]- Grid.TX[ YY ][iy_tay]);
  const double ax = (Grid.XX[0]-Inv.x)/(Grid.XX[0]);
  
  const double xb=Grid.XX[0];
  const double xbsq=xb*xb;
  
  const double ell2a = lerp( ay , Grid.TX[ Gl2 ][iy_tay] , Grid.TX[ Gl2 ][iy2_tay] ) ;
  const double ell2b = accessv(false, true, ix, iy, QL2, false, d0cb, Inv.cb, Inv.y, Grid );
  T -> ell2 = lerp(ax, ell2a, ell2b) ;
  
  const double dell2adx = lerp( ay , Grid.TX[ Gl21 ][iy_tay] , Grid.TX[ Gl21 ][iy2_tay] )*Inv.cb ;
  const double dell2bdx = accessv(false, true, ix, iy, dxQL2, false, d0cb, Inv.cb, Inv.y , Grid );
  T -> dell2dx = lerp(ax, dell2adx, dell2bdx) ;

  const double dell2ady = lerp( ay , Grid.TX[ Gl2dy ][iy_tay] , Grid.TX[ Gl2dy ][iy2_tay] ) ;
  const double dell2bdy = accessv(false, false, ix, iy, QL2 , true, d0cb, Inv.cb, Inv.y , Grid ) ;
  T -> dell2dy = lerp(ax, dell2ady, dell2bdy) ;

  const double dell2adcb = 0.0;
  const double dell2bdcb = accessv(false, true, ix, iy, QL2 , false, d1cb, Inv.cb, Inv.y , Grid );
  T -> dell2dcb = lerp(ax ,dell2adcb ,dell2bdcb ) ;
  
  const double v1b       = accessv(false, true , ix, iy, QG1  , false, d0cb , Inv.cb, Inv.y , Grid );
  const double dv1bdx    = accessv(false, true , ix, iy, dxQG1, false, d0cb , Inv.cb, Inv.y , Grid );
  const double dv1bdy    = accessv(false, false, ix, iy, QG1  , true , d0cb , Inv.cb, Inv.y , Grid );
  const double dv1bdcb   = accessv(false, true , ix, iy, QG1  , false, d1cb , Inv.cb, Inv.y , Grid );
  const double ell4b     = accessv(false, true , ix, iy, QL4  , false, d0cb , Inv.cb, Inv.y , Grid );
  const double dell4bdx  = accessv(false, true , ix, iy, dxQL4, false, d0cb , Inv.cb, Inv.y , Grid );
  const double dell4bdy  = accessv(false, false, ix, iy, QL4  , true , d0cb , Inv.cb, Inv.y , Grid );
  const double dell4bdcb = accessv(false, true , ix, iy, QL4  , false, d1cb , Inv.cb, Inv.y , Grid );
  
  const double ell3a = lerp( ay , Grid.TX[ Gl3 ][iy_tay] , Grid.TX[ Gl3 ][iy2_tay] ) ;
  const double ell3b = ell4b/(2.0*xbsq*ysq) - Inv.y*Inv.cb*ell2b/xb;
  T -> ell3 = lerp( ax , ell3a , ell3b ) ;

  const double ell1b = 4.0/(3.0*xbsq*xbsq)*(v1b - xbsq*ysq*(Inv.cb*Inv.cb-0.25)*ell2b - 1.5*xbsq*xb*Inv.y*Inv.cb*ell3b);
  T -> ell1 = ell1b;

  const double dell3ady = lerp( ay , Grid.TX[ Gl3dy ][iy_tay] , Grid.TX[ Gl3dy ][iy2_tay] ) ;
  const double dell3bdy = (-ell4b/(ysq*Inv.y)
			   + dell4bdy/(2.0*ysq)
			   - Inv.cb*xb*ell2b
			   - Inv.y*xb*Inv.cb*dell2bdy)/(xbsq);
  T -> dell3dy = lerp( ax, dell3ady , dell3bdy ) ; 

  const double dell3adcb = 0.0;
  const double dell3bdcb = (dell4bdcb/(2.0*xb*ysq)
			    -Inv.y*ell2b - Inv.y*Inv.cb*dell2bdcb)/xb;
  T -> dell3dcb = lerp( ax , dell3adcb , dell3bdcb ) ; 

  const double dell3bdx = (-ell4b/(xb*ysq)
			   + dell4bdx/(2.0*ysq)
			   + Inv.y*Inv.cb*ell2b
			   - Inv.y*xb*Inv.cb*dell2bdx)/xbsq;
  T -> dell3dx = dell3bdx;

  T -> dell1dy = 4.0/(3.0*xbsq*xbsq)*(dv1bdy-2.0*xbsq*Inv.y*(Inv.cb*Inv.cb-0.25)*ell2b
				      -xbsq*ysq*(Inv.cb*Inv.cb-0.25)*dell2bdy
				      -1.5*xbsq*xb*Inv.cb*ell3b
				      -1.5*xbsq*xb*Inv.y*Inv.cb*dell3bdy) ;
  
  T -> dell1dcb = 4.0/(3.0*xbsq*xbsq)
    *(dv1bdcb - 2.0*xbsq*ysq*Inv.cb*ell2b-xbsq*ysq*(Inv.cb*Inv.cb-0.25)*dell2bdcb
      -1.5*xbsq*xb*Inv.y*ell3b-1.5*xbsq*xb*Inv.y*Inv.cb*dell3bdcb);

  T -> dell1dx = -4.0*ell1b/xb + 4.0/(3.0*xbsq*xbsq)
    *(dv1bdx - 2.0*xb*ysq*(Inv.cb*Inv.cb-0.25)*ell2b-xbsq*ysq*(Inv.cb*Inv.cb-0.25)*dell2bdx
      -4.5*xbsq*Inv.y*Inv.cb*ell3b - 1.5*xbsq*xb*Inv.y*Inv.cb*dell3bdx);

#ifdef VERBOSE
  printf( "ay %1.15e\n" , ay ) ;
  printf( "ax %1.15e\n" , ax ) ;
  printf( "xb %1.15e\n" , xb ) ;
  printf( "xbsq %1.15e\n" , xbsq ) ;
  printf( "ysq %1.15e %1.15e\n" , Inv.ysq , ysq ) ;
  printf( "ell2a %1.15e\n" , ell2a ) ;
  printf( "ell2b %1.15e\n" , ell2b ) ;
  printf( "dell2 tmps %1.15e %1.15e\n" , dell2ady , dell2bdy ) ;
  printf( "v1b %1.15e\n" , v1b ) ;
  printf( "dv1bdx %1.15e\n" , dv1bdx ) ;
  printf( "dv1bdy %1.15e\n" , dv1bdy ) ;
  printf( "dv1bdcb %1.15e\n" , dv1bdcb ) ;
  printf( "ell4b %1.15e\n" , ell4b ) ;
  printf( "dell4bdx %1.15e\n" , dell4bdx ) ;
  printf( "dell4bdy %1.15e\n" , dell4bdy ) ;
  printf( "dell4bdcb %1.15e\n" , dell4bdcb ) ;
  printf( "ell3a %1.15e\n" , ell3a ) ;
  printf( "ell3b %1.15e\n" , ell3b ) ;
  printf( "T->ell1 %1.15e\n" , T -> ell1 ) ;
  printf( "T->ell2 %1.15e\n" , T -> ell2 ) ;
  printf( "T->ell3 %1.15e\n" , T -> ell3 ) ;
  printf( "T->dell2dx %1.15e\n" , T -> dell2dx ) ;
  printf( "T->dell2dy %1.15e\n" , T -> dell2dy ) ;
  printf( "T->dell2dcb %1.15e\n" , T -> dell2dcb ) ;
  printf( "T->dell3dy %1.15e\n" , T -> dell3dy ) ;
  printf( "T->dell3dcb %1.15e\n" , T -> dell3dcb ) ;
  printf( "T->dell3dx %1.15e\n" , T -> dell3dx ) ;
  printf( "T->dell1dy %1.15e\n" , T -> dell1dy ) ;
  printf( "T->dell1dcb %1.15e\n" , T -> dell1dcb ) ;
  printf( "T->dell1dx %1.15e\n" , T -> dell1dx ) ;
#endif
  
  return 0 ;
}

static void
XMYSWAP_T_contrib( struct ttmps *T ,
		   const struct invariants Inv )
{
  const double ell2s = T -> ell2;
  const double dells2dx  = T -> dell2dx;
  const double dells2dcb = T -> dell2dcb;
  const double dells2dy  = T -> dell2dy;
  const double cd = ( Inv.yorig-Inv.x*Inv.cborig)/Inv.xmy;
  T->ell2     = ell2s;
  T->dell2dx  = dells2dx + (1.0-Inv.cb*Inv.cb)/Inv.xmy*dells2dcb + Inv.cb*dells2dy;
  T->dell2dcb = -(Inv.ysq*cd*dells2dcb + Inv.x*Inv.yorig*Inv.xmy*dells2dy)/Inv.xmysq;
  T->dell2dy  = -(Inv.cborig+Inv.cb*cd)/Inv.xmy*dells2dcb+cd*dells2dy;
  
  const double ell3s     = -T->ell3;
  const double dells3dx  = -T->dell3dx;
  const double dells3dcb = -T->dell3dcb;
  const double dells3dy  = -T->dell3dy;

  T->ell3    =  ell3s - T->ell2; 
  T->dell3dx =  dells3dx + (1.0-Inv.cb*Inv.cb)/Inv.xmy*dells3dcb + Inv.cb*dells3dy - T->dell2dx;
  T->dell3dcb =  -(Inv.ysq*cd*dells3dcb + Inv.x*Inv.yorig*Inv.xmy*dells3dy)/Inv.xmysq - T->dell2dcb;
  T->dell3dy  =  -(Inv.cborig+Inv.cb*cd)/Inv.xmy*dells3dcb+cd*dells3dy - T->dell2dy;
  
  const double ell1s = T->ell1;
  const double dells1dx  = T->dell1dx;
  const double dells1dcb = T->dell1dcb;
  const double dells1dy  = T->dell1dy;
  T->ell1     = ell1s - T->ell2 -2.0*T->ell3;
  T->dell1dx  = dells1dx + (1.0-Inv.cb*Inv.cb)/Inv.xmy*dells1dcb
    + Inv.cb*dells1dy - T->dell2dx - 2.0*T->dell3dx;
  T->dell1dcb = -(Inv.ysq*cd*dells1dcb + Inv.x*Inv.yorig*Inv.xmy*dells1dy)/Inv.xmysq
    - T->dell2dcb -2.0*T->dell3dcb;
  T->dell1dy  = -(Inv.cborig+Inv.cb*cd)/Inv.xmy*dells1dcb+cd*dells1dy
    - T->dell2dy -2.0*T->dell3dy;
}

// set dyv[alf][bet][dta]= \partial^(y)_\beta <(epsilon_alf epsilon_dta - 1/4 delta_{alf dta}) I>_epsilon
// and dxv[alf][bet][dta]= \partial^(x)_\beta <(epsilon_alf epsilon_dta - 1/4 delta_{alf dta}) I>_epsilon
int
chnr_dT( const double xv[4] ,
	 const double yv[4] ,
	 const struct invariants Inv ,
	 const struct Grid_coeffs Grid ,
	 const struct AVX_precomps *PC ,
	 double dxv[4][4][4] ,
	 double dyv[4][4][4] )
{
  struct ttmps T ;
  const double ysq = Inv.flag2 ? Inv.xmysq : Inv.ysq ;
  
  if( Inv.x < Grid.XX[0]) {     
    if( TAYLORX_contribution( &T , Inv , Grid ) == 1 ) {
      return 1 ;
    }
  } else if( Inv.x > Grid.XX[ Grid.nstpx -1 ] ) {
    return 1 ;
  } else {
    double f[4] KQED_ALIGN ;
    extractff2( QL2 , d0cb , Inv , Grid , PC , f ) ;
    T.dell2dx  = f[1] ;
    T.dell2dy  = f[2] ;
    T.ell2     = f[3] ;
    T.dell2dcb = f[0] ;

    extractff2( QL4 , d0cb , Inv , Grid , PC , f ) ;
    T.dell4dx  = f[1] ;
    T.dell4dy  = f[2] ;
    T.ell4     = f[3] ;
    T.dell4dcb = f[0] ; 
    
    extractff2( QG1 , d0cb , Inv , Grid , PC , f ) ;
    T.dv1dx  = f[1] ;
    T.dv1dy  = f[2] ;
    T.v1     = f[3] ;
    T.dv1dcb = f[0] ;
    
    T.ell3 = T.ell4/(2.0*Inv.xsq*ysq) - Inv.y*Inv.cb*T.ell2/Inv.x;
    T.ell1 = 4.0/(3.0*Inv.xsq*Inv.xsq)*
      (T.v1 - Inv.xsq*ysq*(Inv.cb*Inv.cb-0.25)*T.ell2
       - 1.5*Inv.xsq*Inv.x*Inv.y*Inv.cb*T.ell3);
        
    T.dell3dy  = (-T.ell4/(ysq*Inv.y) + T.dell4dy/(2.0*ysq)
		  - Inv.cb*Inv.x*T.ell2 -
		  Inv.y*Inv.x*Inv.cb*T.dell2dy)/(Inv.xsq);
    
    T.dell3dcb = (T.dell4dcb/(2.0*Inv.x*ysq) -
		  Inv.y*T.ell2 - Inv.y*Inv.cb*T.dell2dcb)/Inv.x;
    
    T.dell3dx  = (-T.ell4/(Inv.x*ysq)
		  + T.dell4dx/(2.0*ysq)
		  + Inv.y*Inv.cb*T.ell2
		  - Inv.y*Inv.x*Inv.cb*T.dell2dx)/Inv.xsq;
    
    T.dell1dy  = 4.0/(3.0*Inv.xsq*Inv.xsq)*(T.dv1dy-2.0*Inv.xsq*Inv.y*(Inv.cb*Inv.cb-0.25)*T.ell2
				    -Inv.xsq*ysq*(Inv.cb*Inv.cb-0.25)*T.dell2dy
				    -1.5*Inv.xsq*Inv.x*Inv.cb*T.ell3-1.5*Inv.xsq*Inv.x*Inv.y*Inv.cb*T.dell3dy);
    
    T.dell1dcb = 4.0/(3.0*Inv.xsq*Inv.xsq)*(T.dv1dcb
				    - 2.0*Inv.xsq*ysq*Inv.cb*T.ell2-Inv.xsq*ysq*(Inv.cb*Inv.cb-0.25)*T.dell2dcb
				    -1.5*Inv.xsq*Inv.x*Inv.y*T.ell3-1.5*Inv.xsq*Inv.x*Inv.y*Inv.cb*T.dell3dcb);
    
    T.dell1dx  = -4.0*T.ell1/Inv.x
      + 4.0/(3.0*Inv.xsq*Inv.xsq)*(T.dv1dx- 2.0*Inv.x*ysq*(Inv.cb*Inv.cb-0.25)*T.ell2
			   -Inv.xsq*ysq*(Inv.cb*Inv.cb-0.25)*T.dell2dx
			   -4.5*Inv.xsq*Inv.y*Inv.cb*T.ell3
				   - 1.5*Inv.xsq*Inv.x*Inv.y*Inv.cb*T.dell3dx);
  }

  if( Inv.flag2 ) {
    XMYSWAP_T_contrib( &T , Inv ) ;
  } 
   
  int alf, bet, dta ;
  for(bet=0;bet<4;bet++)  {  
    const double c1 = yv[bet]/Inv.yorig;
    const double c3 = xv[bet]/Inv.x;
    const double c2 = (c3-Inv.cborig*c1)/Inv.yorig;
    const double c4 = (c1-Inv.cborig*c3)/Inv.x;
    for(alf=0;alf<4;alf++)  {
      const int d_ab = ( alf==bet ) ;
      for(dta=0;dta<4;dta++) {
	const int d_bd = ( bet==dta ) ;
	const int d_ad = ( alf==dta ) ;
	const double t1 = d_ab*yv[dta]+d_bd*yv[alf]-0.5*d_ad*yv[bet];
	const double t2 = d_ab*xv[dta]+d_bd*xv[alf]-0.5*d_ad*xv[bet];
	const double t3 = xv[alf]*xv[dta]-d_ad*Inv.xsq/4;
	const double t4 = yv[alf]*yv[dta]-d_ad*Inv.ysq/4;
	const double t5 = xv[alf]*yv[dta]+yv[alf]*xv[dta]-0.5*Inv.xdoty*d_ad;
	dyv[alf][bet][dta] =
	  + t1*T.ell2
	  + t2*T.ell3
	  + t3*(c1*T.dell1dy + c2*T.dell1dcb)
	  + t4*(c1*T.dell2dy + c2*T.dell2dcb)
	  + t5*(c1*T.dell3dy + c2*T.dell3dcb)
	  ;
	dxv[alf][bet][dta] =
	  + t2*T.ell1
	  + t1*T.ell3
	  + t3*(c3*T.dell1dx + c4*T.dell1dcb)
	  + t4*(c3*T.dell2dx + c4*T.dell2dcb)
	  + t5*(c3*T.dell3dx + c4*T.dell3dcb)
	  ;
      }
    }
  }
  
  return 0 ;
}
