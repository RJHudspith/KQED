/**
   @file Tabd.c
   @brief contains the functions Tabd_xeq0 and Tabd_yeq0

   @warning these differ by a factor of 2 from the old code as they contain the implicit SYMG normalisation
 */
#include "KQED.h"

#include "chnr_dV.h"   // lerp
#include "getff-new.h" // SCALPROD, bsrch

// delta function
static const double dlta[4][4] = { { 1 , 0 , 0 , 0 } ,
				   { 0 , 1 , 0 , 0 } ,
				   { 0 , 0 , 1 , 0 } ,
				   { 0 , 0 , 0 , 1 } } ;

// unroll some inner loops in the hope that it speeds things up
#define UNROLL_BETDTA

// Taylor expanded versions for x == 0
int
Tabd_xeq0( const double yv[4] ,
	   const struct Grid_coeffs Grid ,
	   double tI[4][4][4] ,
	   double tII[4][4][4] ,
	   double tIII[4][4][4] )
{ 
  const double ysq = SCALPROD(yv,yv);
  const double y = sqrt( fabs( ysq ) );
  
  const int iy_tay = bsrch( Grid.TX[ YY ] , y , 0 , Grid.NY_tay ) ;

  if( iy_tay == Grid.NY_tay-1 ) {
    #ifdef verbose
    fprintf( stderr , "# Warning in Tabd_xeq0: y is higher than upper edge.\n");
    fprintf( stderr , "Setting kernel to ZERO (%f > %f)\n" ,
	     y , Grid.TX[ YY ][ Grid.NY_tay-1] ) ;
    #endif
    // need to return here otherwise we will segfault
    return 1 ;
  }
  // can set this after we checked the above error
  const int iy2_tay = iy_tay+1 ;
  
  const double ay = (Grid.TX[ YY ][iy2_tay]-y)/( Grid.TX[ YY ][iy2_tay]-Grid.TX[ YY ][iy_tay]);
  
  const double fx       = 0.5*lerp( ay , Grid.TX[ G0dx ][iy_tay] , Grid.TX[ G0dx ][iy2_tay] ) ; 
  const double fy       = 0.5*lerp( ay , Grid.TX[ G0dy ][iy_tay] , Grid.TX[ G0dy ][iy2_tay] ) ;
  const double yhat[4]  = { yv[0]/y , yv[1]/y , yv[2]/y , yv[3]/y } ;
  
  const double dg0dy[4] = { yhat[0]*fy/2 , yhat[1]*fy/2 , yhat[2]*fy/2 , yhat[3]*fy/2 } ; 
  const double dg0dx[4] = { yhat[0]*fx/2 , yhat[1]*fx/2 , yhat[2]*fx/2 , yhat[3]*fx/2 } ;
  
  const double ell2     = lerp( ay , Grid.TX[ Gl2 ][iy_tay] , Grid.TX[ Gl2 ][iy2_tay] ) ; 
  const double dell2adx = lerp( ay , Grid.TX[ Gl21 ][iy_tay] , Grid.TX[ Gl21 ][iy2_tay] ) ;
  const double dell2ady = lerp( ay , Grid.TX[ Gl2dy ][iy_tay] , Grid.TX[ Gl2dy ][iy2_tay] ) ;
  const double ell3a    = lerp( ay , Grid.TX[ Gl3 ][iy_tay] , Grid.TX[ Gl3 ][iy2_tay] ) ;
  
  const double dg2dx    = lerp( ay , Grid.TX[ G21 ][iy_tay] , Grid.TX[ G21 ][iy2_tay] ) ;
  const double ddg2dxdy = lerp( ay , Grid.TX[ G21dy ][iy_tay] , Grid.TX[ G21dy ][iy2_tay] ) ;
  const double dg1dx    = (lerp( ay , Grid.TX[ G31A ][iy_tay] , Grid.TX[ G31A ][iy2_tay] )
			   -(2./3)*( lerp( ay , Grid.TX[ G31B ][iy_tay] , Grid.TX[ G31B ][iy2_tay] ) )
			   -y* ( lerp( ay , Grid.TX[ G22A ][iy_tay] , Grid.TX[ G22A ][iy2_tay] ) ) ) ;
  const double dg1dy    = (lerp( ay , Grid.TX[ G3Ady ][iy_tay] , Grid.TX[ G3Ady ][iy2_tay] )
			   -lerp( ay , Grid.TX[ G3Bdy ][iy_tay] , Grid.TX[ G3Bdy ][iy2_tay] ) ) ;
  const double phi1     = ( lerp( ay , Grid.TX[ G22A ][iy_tay] , Grid.TX[ G22A ][iy2_tay] )
			    - lerp( ay , Grid.TX[ G22B ][iy_tay] , Grid.TX[ G22B ][iy2_tay] )/5. ) ;
  const double phi2     = ( lerp( ay , Grid.TX[ G22B ][iy_tay] , Grid.TX[ G22B ][iy2_tay] ) )*1.2 ;
  
  // various precomputations to make the loops faster
  const double ysq_4 = ysq/4. ;
  register const double E = (ell3a + ell2) , D2 = (dell2adx + dell2ady) ,
    D1 = (dg1dy+dg1dx) , D12 = (dg1dx+dg2dx) , D3 = (y*ddg2dxdy-dg2dx) ;
  const double yvE[4]       = { yv[0]*E , yv[1]*E , yv[2]*E , yv[3]*E } ;
  const double yvell3a[4]   = { yv[0]*ell3a , yv[1]*ell3a , yv[2]*ell3a , yv[3]*ell3a } ;
  const double yhatD1[4]    = { yhat[0]*D1 , yhat[1]*D1 , yhat[2]*D1 , yhat[3]*D1 } ;
  const double yhatD2[4]    = { yhat[0]*D2 , yhat[1]*D2 , yhat[2]*D2 , yhat[3]*D2 } ;
  const double yhatD3[4]    = { yhat[0]*D3 , yhat[1]*D3 , yhat[2]*D3 , yhat[3]*D3 } ;
  const double yhatdg2dx[4] = { yhat[0]*dg2dx , yhat[1]*dg2dx , yhat[2]*dg2dx , yhat[3]*dg2dx } ;
  
  // point out the structures
  register double *PtIII = (double*)tIII ;
  register double *PtII  = (double*)tII ;
  register double *PtI   = (double*)tI ;
  
  int alf , dta , bet ;
  for(alf=0;alf<4;alf++) {
    register const double yalf2dx = yhat[alf]*(dell2adx) ;
    register const double tIIsum  = dg0dx[alf]-0.5*yvell3a[alf]-ysq_4*yalf2dx;
    register const double yD12    = yhat[alf]*D12 ;
    for(bet=0;bet<4;bet++) {
      register const double dgSUM   = (dg0dx[bet]+dg0dy[bet]) ; 
      register const double tIsum   = 2*(dlta[alf][bet]*phi1+yhat[alf]*yhat[bet]*phi2) ;
      register const double tIIIsum = dgSUM-0.5*yvE[bet]-ysq_4*yhatD2[bet] ;
      for(dta=0;dta<4;dta++) {
	// (d/dxbeta+d/dybeta) T_{alpha delta}
	*PtIII = (yv[alf]*yv[dta]*yhatD2[bet]
		  +dlta[alf][bet]*yvE[dta]
		  +dlta[bet][dta]*yvE[alf]
		  +dlta[alf][dta]*tIIIsum) ; PtIII++ ;
	// d/dxalpha T_{beta delta}
	*PtII = (yv[bet]*yv[dta]*yalf2dx+dlta[bet][dta]*(tIIsum) 
		 +(dlta[bet][alf]*yvell3a[dta]+dlta[alf][dta]*yvell3a[bet])) ; PtII++ ;
	// first the d/dxalpha d/dybeta V_delta terms:  
        *PtI = (dlta[bet][dta]*yD12+yhatdg2dx[dta]*dlta[alf][bet]
		+yhatD3[bet]*yhat[alf]*yhat[dta]
		+yhatD1[bet]*dlta[alf][dta]+yv[dta]*tIsum) ; PtI++ ;
      }
    }
  }
    
  return 0 ;
}

// Taylor expansions around y = 0
int
Tabd_yeq0( const double xv[4] ,
	   const struct Grid_coeffs Grid ,
	   double tI[4][4][4] ,
	   double tII[4][4][4] ,
	   double tIII[4][4][4] )
{  
  const double xsq = SCALPROD(xv,xv);
  const double x = sqrt( fabs(xsq) ); 

  const int ix_tay = bsrch( Grid.XX , x , 0 , Grid.nstpx ) ;
  
  if( x > Grid.XX[ Grid.nstpx - 1 ] ) {
    #ifdef verbose
    fprintf( stderr , "# Warning in Tabd_yeq0: x is higher than upper edge.\n");
    fprintf( stderr , "Setting kernel to ZERO (%f > %f)\n" ,
	     x , Grid.XX[ Grid.nstpx-1] ) ;
    #endif
    // need to return here otherwise we will segfault
    return 1 ;
  }
  // we can set this now as we aren't at the upper edge
  const int ix2_tay = ix_tay+1 ;
  
  const double ax = ( Grid.XX[ix2_tay]-x)/( Grid.XX[ix2_tay]-Grid.XX[ix_tay]);
  const double xa = Grid.XX[ix_tay];
  const double xb = Grid.XX[ix2_tay];
  const double fx = 0.5*lerp( ax , Grid.TY[ alpha0dx_0p ][ ix_tay ], Grid.TY[ alpha0dx_0p ][ix2_tay] ) ;
  const double fy = ( lerp( ax , Grid.TY[ alpha0_1p ][ix_tay] , Grid.TY[ alpha0_1p ][ix2_tay] ) );

  // set by hand the array variables
  const double xhat[4] = { xv[0]/x , xv[1]/x , xv[2]/x , xv[3]/x } ;
  const double dg0dy[4] = { xhat[0]*fy/2 , xhat[1]*fy/2 , xhat[2]*fy/2 , xhat[3]*fy/2 } ;
  const double dg0dx[4] = { xhat[0]*fx/2 , xhat[1]*fx/2 , xhat[2]*fx/2 , xhat[3]*fx/2 } ;

  const double xa2 = xa*xa ;
  const double xa3 = xa2*xa ;
  const double xa4 = xa3*xa ;
  const double xb2 = xb*xb ;
  const double xb3 = xb2*xb ;
  const double xb4 = xb3*xb ;

  const double ell1       = 4*( lerp( ax , Grid.TY[alpha1_0p][ix_tay]/(xa4) ,
				      Grid.TY[alpha1_0p][ix2_tay]/(xb4) ) )/3. ;
  const double ell3       = lerp( ax , Grid.TY[ beta4_1p ][ix_tay]/(xa2) ,
				  Grid.TY[ beta4_1p ][ix2_tay]/(xb2) ) ;
  
  const double dell1dx    = 4*lerp( ax , (Grid.TY[alpha1dx_0p][ix_tay]
					  -4.*Grid.TY[alpha1_0p][ix_tay]/xa)/xa4 ,
				    (Grid.TY[alpha1dx_0p][ix2_tay]
				     -4.*Grid.TY[alpha1_0p][ix2_tay]/xb)/xb4 )/3. ;
  
  const double dell1dycb  = 2*lerp( ax , (4.*Grid.TY[alpha1_1p][ix_tay]/(3*xa)
					  - Grid.TY[beta4_1p][ix_tay])/xa3 ,
				    (4.*Grid.TY[alpha1_1p][ix2_tay]/(3*xb)
				     - Grid.TY[beta4_1p][ix2_tay])/xb3 ) ;
  
  const double dg1dx      = lerp( ax , Grid.TY[ alpha3dx_0p ][ix_tay] , Grid.TY[ alpha3dx_0p ][ix2_tay] ) ;
  const double ddg1dxdx   = lerp( ax , Grid.TY[ alpha3dxdx_0p ][ix_tay] , Grid.TY[ alpha3dxdx_0p ][ix2_tay] ) ;
  const double dg2dx      = 2*lerp( ax , Grid.TY[ beta2dx_1p ][ix_tay] , Grid.TY[ beta2dx_1p ][ix2_tay] ) ;
  
  const double dg1dycb    = 2*lerp( ax , Grid.TY[ alpha3_1p ][ix_tay]
				    - Grid.TY[beta2_1p ][ix_tay]/xa ,
				    Grid.TY[ alpha3_1p ][ix2_tay]
				    - Grid.TY[beta2_1p ][ix2_tay]/xb ) ;
  
  const double ddg1dxdycb = 2*lerp( ax , Grid.TY[beta2_1p][ix_tay]/(xa2)
				    - Grid.TY[beta2dx_1p][ix_tay]/xa
				    + Grid.TY[alpha3dx_1p][ix_tay] ,
				    Grid.TY[beta2_1p][ix2_tay]/(xb2)
				    - Grid.TY[beta2dx_1p][ix2_tay]/xb
				    + Grid.TY[alpha3dx_1p][ix2_tay] ) ;

  // precomputations
  register const double xsq_4 = xsq/4. , E = ell1+ell3 , D1 = dell1dx + dell1dycb
    , D2 = (ddg1dxdx+ddg1dxdycb) , D3 = (dg1dycb-dg2dx) , D4 = (dg1dx+dg1dycb) ;
  const double xvE[4]    = { xv[0]*E , xv[1]*E , xv[2]*E , xv[3]*E } ;
  const double xvell1[4] = { xv[0]*ell1 , xv[1]*ell1 , xv[2]*ell1 , xv[3]*ell1 } ;

  // point out the structures
  register double *PtIII = (double*)tIII ;
  register double *PtII  = (double*)tII ;
  register double *PtI   = (double*)tI ;

  // I found that unrolling the bet/det loop and performing some CSE helped a little
  // not really a great deal but this is the hottest part of this code so every little helps
#ifndef UNROLL_BETDTA
  int alf , dta , bet ;
  for(alf=0;alf<4;alf++) {
    const register double xhd1dx = xhat[alf]*(dell1dx) ;
    register const double tIIsum = (dg0dx[alf]-0.5*xvell1[alf] ); 
    for(bet=0;bet<4;bet++) {
      register const double Dg0SUM = (dg0dx[bet]+dg0dy[bet]-0.5*xvE[bet]) ;
      for(dta=0;dta<4;dta++) {
	register const double xad = (xv[alf]*xv[dta]-xsq_4*dlta[alf][dta]);
	*PtIII = ((xad*xhat[bet]*D1 + dlta[alf][bet]*xvE[dta]+dlta[bet][dta]*xvE[alf]+Dg0SUM*dlta[alf][dta]) ); PtIII++ ;
	*PtII  = ((xv[bet]*xv[dta]-xsq_4*dlta[bet][dta])*xhd1dx+dlta[bet][alf]*xvell1[dta]+dlta[alf][dta]*xvell1[bet]+dlta[bet][dta]*(tIIsum)); PtII++ ;
	*PtI   = ((dlta[alf][dta]*xhat[bet]+dlta[bet][dta]*xhat[alf]+dlta[alf][bet]*xhat[dta]-xhat[alf]*xhat[bet]*xhat[dta])*D4+xhat[alf]*( xhat[bet]*xv[dta]*D2 - dlta[bet][dta]*D3 ));PtI++ ;	
      }
    }
  }
#else
  register double Dg0SUM , xad , xhd1dx , tIIsum ;
  int alf ;
  for( alf = 0 ; alf < 4 ; alf++ ) {
    xhd1dx = xhat[alf]*(dell1dx) ;
    tIIsum = dg0dx[alf]-0.5*xvell1[alf] ; 
    // bet = 0
    //  dta = 0
    Dg0SUM = (dg0dx[0]+dg0dy[0]-0.5*xvE[0]);
    xad    = (xv[alf]*xv[0]-xsq_4*dlta[alf][0]);
    *PtIII = (xad*xhat[0]*D1 + dlta[alf][0]*xvE[0]+xvE[alf]+Dg0SUM*dlta[alf][0] ); PtIII++ ;
    *PtII  = ((xv[0]*xv[0]-xsq_4)*xhd1dx+dlta[0][alf]*xvell1[0]+dlta[alf][0]*xvell1[0]+tIIsum) ; PtII++ ;
    *PtI   = ((dlta[alf][0]*xhat[0]+xhat[alf]+dlta[alf][0]*xhat[0]-xhat[alf]*xhat[0]*xhat[0])*D4+xhat[alf]*(xhat[0]*xv[0]*D2-D3));PtI++ ;
    //  dta = 1
    xad    = (xv[alf]*xv[1]-xsq_4*dlta[alf][1]);
    *PtIII = (xad*xhat[0]*D1+dlta[alf][0]*xvE[1]+Dg0SUM*dlta[alf][1] ); PtIII++ ;
    *PtII  = ((xv[0]*xv[1])*xhd1dx+dlta[0][alf]*xvell1[1]+dlta[alf][1]*xvell1[0] ); PtII++ ;
    *PtI   = ((dlta[alf][1]*xhat[0]+dlta[alf][0]*xhat[1]-xhat[alf]*xhat[0]*xhat[1])*D4+xhat[alf]*(xhat[0]*xv[1]*D2));PtI++ ;
    //  dta = 2
    xad    = (xv[alf]*xv[2]-xsq_4*dlta[alf][2]);
    *PtIII = (xad*xhat[0]*D1 + dlta[alf][0]*xvE[2]+Dg0SUM*dlta[alf][2] ); PtIII++ ;
    *PtII  = ((xv[0]*xv[2]-xsq_4*dlta[0][2])*xhd1dx+dlta[0][alf]*xvell1[2]+dlta[alf][2]*xvell1[0] ) ; PtII++ ;
    *PtI   = ((dlta[alf][2]*xhat[0]+dlta[alf][0]*xhat[2]-xhat[alf]*xhat[0]*xhat[2])*D4+xhat[alf]*(xhat[0]*xv[2]*D2)); PtI++ ;	
    //  dta = 3
    xad    = (xv[alf]*xv[3]-xsq_4*dlta[alf][3]);
    *PtIII = (xad*xhat[0]*D1 + dlta[alf][0]*xvE[3]+Dg0SUM*dlta[alf][3]) ; PtIII++ ;
    *PtII  = ((xv[0]*xv[3])*xhd1dx+dlta[0][alf]*xvell1[3]+dlta[alf][3]*xvell1[0]) ; PtII++ ;
    *PtI   = ((dlta[alf][3]*xhat[0]+dlta[alf][0]*xhat[3]-xhat[alf]*xhat[0]*xhat[3])*D4+xhat[alf]*(xhat[0]*xv[3]*D2)) ;PtI++ ;
    // bet = 1
    //  dta = 0
    Dg0SUM = (dg0dx[1]+dg0dy[1]-0.5*xvE[1]) ;      
    xad    = (xv[alf]*xv[0]-xsq_4*dlta[alf][0]);
    *PtIII = (xad*xhat[1]*D1 + dlta[alf][1]*xvE[0]+Dg0SUM*dlta[alf][0]) ; PtIII++ ;
    *PtII  = ((xv[1]*xv[0])*xhd1dx+dlta[1][alf]*xvell1[0]+dlta[alf][0]*xvell1[1]) ; PtII++ ;
    *PtI   = ((dlta[alf][0]*xhat[1]+dlta[alf][1]*xhat[0]-xhat[alf]*xhat[1]*xhat[0])*D4+xhat[alf]*(xhat[1]*xv[0]*D2)) ;PtI++ ;
    //  dta = 1
    xad    = (xv[alf]*xv[1]-xsq_4*dlta[alf][1]);
    *PtIII = (xad*xhat[1]*D1 + dlta[alf][1]*xvE[1]+xvE[alf]+Dg0SUM*dlta[alf][1] ); PtIII++ ;
    *PtII  = ((xv[1]*xv[1]-xsq_4)*xhd1dx+dlta[1][alf]*xvell1[1]+dlta[alf][1]*xvell1[1]+tIIsum ); PtII++ ;
    *PtI   = ((dlta[alf][1]*xhat[1]+xhat[alf]+dlta[alf][1]*xhat[1]-xhat[alf]*xhat[1]*xhat[1])*D4+xhat[alf]*(xhat[1]*xv[1]*D2-D3));PtI++ ;
    //  dta = 2
    xad    = (xv[alf]*xv[2]-xsq_4*dlta[alf][2]);
    *PtIII = (xad*xhat[1]*D1 + dlta[alf][1]*xvE[2]+Dg0SUM*dlta[alf][2]) ; PtIII++ ;
    *PtII  = ((xv[1]*xv[2])*xhd1dx+dlta[1][alf]*xvell1[2]+dlta[alf][2]*xvell1[1] ); PtII++ ;
    *PtI   = ((dlta[alf][2]*xhat[1]+dlta[alf][1]*xhat[2]-xhat[alf]*xhat[1]*xhat[2])*D4+xhat[alf]*(xhat[1]*xv[2]*D2)) ;PtI++ ;	
    //  dta = 3
    xad    = (xv[alf]*xv[3]-xsq_4*dlta[alf][3]);
    *PtIII = (xad*xhat[1]*D1 + dlta[alf][1]*xvE[3]+Dg0SUM*dlta[alf][3] ); PtIII++ ;
    *PtII  = ((xv[1]*xv[3])*xhd1dx+dlta[1][alf]*xvell1[3]+dlta[alf][3]*xvell1[1] ); PtII++ ;
    *PtI   = ((dlta[alf][3]*xhat[1]+dlta[alf][1]*xhat[3]-xhat[alf]*xhat[1]*xhat[3])*D4+xhat[alf]*(xhat[1]*xv[3]*D2));PtI++ ;
    // bet = 2
    //  dta = 0
    Dg0SUM = (dg0dx[2]+dg0dy[2]-0.5*xvE[2]);      
    xad    = (xv[alf]*xv[0]-xsq_4*dlta[alf][0]);
    *PtIII = (xad*xhat[2]*D1 + dlta[alf][2]*xvE[0]+Dg0SUM*dlta[alf][0] ); PtIII++ ;
    *PtII  = ((xv[2]*xv[0])*xhd1dx+dlta[2][alf]*xvell1[0]+dlta[alf][0]*xvell1[2] ); PtII++ ;
    *PtI   = ((dlta[alf][0]*xhat[2]+dlta[alf][2]*xhat[0]-xhat[alf]*xhat[2]*xhat[0])*D4+xhat[alf]*(xhat[2]*xv[0]*D2));PtI++ ;
    //  dta = 1
    xad    = (xv[alf]*xv[1]-xsq_4*dlta[alf][1]);
    *PtIII = (xad*xhat[2]*D1 + dlta[alf][2]*xvE[1]+Dg0SUM*dlta[alf][1] ); PtIII++ ;
    *PtII  = ((xv[2]*xv[1])*xhd1dx+dlta[2][alf]*xvell1[1]+dlta[alf][1]*xvell1[2] ) ; PtII++ ;
    *PtI   = ((dlta[alf][1]*xhat[2]+dlta[alf][2]*xhat[1]-xhat[alf]*xhat[2]*xhat[1])*D4+xhat[alf]*(xhat[2]*xv[1]*D2));PtI++ ;
    //  dta = 2
    xad    = (xv[alf]*xv[2]-xsq_4*dlta[alf][2]);
    *PtIII = (xad*xhat[2]*D1 + dlta[alf][2]*xvE[2]+xvE[alf]+Dg0SUM*dlta[alf][2] ); PtIII++ ;
    *PtII  = ((xv[2]*xv[2]-xsq_4)*xhd1dx+dlta[2][alf]*xvell1[2]+dlta[alf][2]*xvell1[2]+tIIsum ); PtII++ ;
    *PtI   = ((dlta[alf][2]*xhat[2]+xhat[alf]+dlta[alf][2]*xhat[2]-xhat[alf]*xhat[2]*xhat[2])*D4+xhat[alf]*(xhat[2]*xv[2]*D2-D3));PtI++ ;
    //  dta = 3
    xad    = (xv[alf]*xv[3]-xsq_4*dlta[alf][3]);
    *PtIII = (xad*xhat[2]*D1 + dlta[alf][2]*xvE[3]+Dg0SUM*dlta[alf][3] ); PtIII++ ;
    *PtII  = ((xv[2]*xv[3])*xhd1dx+dlta[2][alf]*xvell1[3]+dlta[alf][3]*xvell1[2] ); PtII++ ;
    *PtI   = ((dlta[alf][3]*xhat[2]+dlta[alf][2]*xhat[3]-xhat[alf]*xhat[2]*xhat[3])*D4+xhat[alf]*(xhat[2]*xv[3]*D2));PtI++ ;
    // bet = 3
    //  dta = 0
    Dg0SUM = (dg0dx[3]+dg0dy[3]-0.5*xvE[3]);      
    xad    = (xv[alf]*xv[0]-xsq_4*dlta[alf][0]);
    *PtIII = (xad*xhat[3]*D1 + dlta[alf][3]*xvE[0]+Dg0SUM*dlta[alf][0] ); PtIII++ ;
    *PtII  = ((xv[3]*xv[0])*xhd1dx+dlta[3][alf]*xvell1[0]+dlta[alf][0]*xvell1[3] ) ; PtII++ ;
    *PtI   = ((dlta[alf][0]*xhat[3]+dlta[alf][3]*xhat[0]-xhat[alf]*xhat[3]*xhat[0])*D4+xhat[alf]*(xhat[3]*xv[0]*D2)); PtI++ ;
    //  dta = 1
    xad    = (xv[alf]*xv[1]-xsq_4*dlta[alf][1]);
    *PtIII = (xad*xhat[3]*D1 + dlta[alf][3]*xvE[1]+Dg0SUM*dlta[alf][1] ); PtIII++ ;
    *PtII  = ((xv[3]*xv[1])*xhd1dx+dlta[3][alf]*xvell1[1]+dlta[alf][1]*xvell1[3] ); PtII++ ;
    *PtI   = ((dlta[alf][1]*xhat[3]+dlta[alf][3]*xhat[1]-xhat[alf]*xhat[3]*xhat[1])*D4+xhat[alf]*(xhat[3]*xv[1]*D2));PtI++ ;
    //  dta = 2
    xad    = (xv[alf]*xv[2]-xsq_4*dlta[alf][2]);
    *PtIII = (xad*xhat[3]*D1 + dlta[alf][3]*xvE[2]+Dg0SUM*dlta[alf][2]); PtIII++ ;
    *PtII  = ((xv[3]*xv[2])*xhd1dx+dlta[3][alf]*xvell1[2]+dlta[alf][2]*xvell1[3] ); PtII++ ;
    *PtI   = ((dlta[alf][2]*xhat[3]+dlta[alf][3]*xhat[2]-xhat[alf]*xhat[3]*xhat[2])*D4+xhat[alf]*(xhat[3]*xv[2]*D2));PtI++ ;	
    //  dta = 3
    xad    = (xv[alf]*xv[3]-xsq_4*dlta[alf][3]);
    *PtIII = (xad*xhat[3]*D1 + dlta[alf][3]*xvE[3]+xvE[alf]+Dg0SUM*dlta[alf][3] ); PtIII++ ;
    *PtII  = ((xv[3]*xv[3]-xsq_4)*xhd1dx+dlta[3][alf]*xvell1[3]+dlta[alf][3]*xvell1[3]+tIIsum ); PtII++ ;
    *PtI   = ((dlta[alf][3]*xhat[3]+xhat[alf]+dlta[alf][3]*xhat[3]-xhat[alf]*xhat[3]*xhat[3])*D4+xhat[alf]*(xhat[3]*xv[3]*D2-D3)); PtI++ ;
  }
#endif
  
  return 0 ;
}
