/**
   @file Tabd.c
   @brief contains the functions Tabd_xeq0 and Tabd_yeq0

   @warning these differ by a factor of 2 from the old code as they contain the implicit SYMG normalisation
 */
#include "KQED.h"

#include "chnr_dV.h"   // lerp
#include "getff-new.h" // SCALPROD, bsrch

// delta function
__device__
static const double dlta[4][4] = { { 1 , 0 , 0 , 0 } ,
				   { 0 , 1 , 0 , 0 } ,
				   { 0 , 0 , 1 , 0 } ,
				   { 0 , 0 , 0 , 1 } } ;

// Taylor expanded versions for x == 0
__device__ KQED_PRIVATE
int
Tabd_xeq0( const double yv[4] ,
	   const struct Grid_coeffs Grid ,
	   double tI[4][4][4] ,
	   double tII[4][4][4] ,
	   double tIII[4][4][4] )
{ 
  const double ysq = SCALPROD(yv,yv);
  const double y = sqrt( fabs( ysq ) );
  
  const int iy_tay = find_ind( getTX(&Grid,YY) , y , 0 , Grid.NY_tay ) ;

  if( iy_tay == Grid.NY_tay-1 ) {
    // need to return here otherwise we will segfault
    return 1 ;
  }
  // can set this after we checked the above error
  const int iy2_tay = iy_tay+1 ;
  
  const double ay = (getTX(&Grid,YY)[iy2_tay]-y)/( getTX(&Grid,YY)[iy2_tay]-getTX(&Grid,YY)[iy_tay]);
  
  const double fx       = 0.5*lerp( ay , getTX(&Grid,G0dx)[iy_tay] , getTX(&Grid,G0dx)[iy2_tay] ) ; 
  const double fy       = 0.5*lerp( ay , getTX(&Grid,G0dy)[iy_tay] , getTX(&Grid,G0dy)[iy2_tay] ) ;
  const double yhat[4]  = { yv[0]/y , yv[1]/y , yv[2]/y , yv[3]/y } ;
  
  const double dg0dy[4] = { yhat[0]*fy/2 , yhat[1]*fy/2 , yhat[2]*fy/2 , yhat[3]*fy/2 } ; 
  const double dg0dx[4] = { yhat[0]*fx/2 , yhat[1]*fx/2 , yhat[2]*fx/2 , yhat[3]*fx/2 } ;
  
  const double ell2     = lerp( ay , getTX(&Grid,Gl2)[iy_tay] , getTX(&Grid,Gl2)[iy2_tay] ) ; 
  const double dell2adx = lerp( ay , getTX(&Grid,Gl21)[iy_tay] , getTX(&Grid,Gl21)[iy2_tay] ) ;
  const double dell2ady = lerp( ay , getTX(&Grid,Gl2dy)[iy_tay] , getTX(&Grid,Gl2dy)[iy2_tay] ) ;
  const double ell3a    = lerp( ay , getTX(&Grid,Gl3)[iy_tay] , getTX(&Grid,Gl3)[iy2_tay] ) ;
  
  const double dg2dx    = lerp( ay , getTX(&Grid,G21)[iy_tay] , getTX(&Grid,G21)[iy2_tay] ) ;
  const double ddg2dxdy = lerp( ay , getTX(&Grid,G21dy)[iy_tay] , getTX(&Grid,G21dy)[iy2_tay] ) ;
  const double dg1dx    = (lerp( ay , getTX(&Grid,G31A)[iy_tay] , getTX(&Grid,G31A)[iy2_tay] )
			   -(2./3)*( lerp( ay , getTX(&Grid,G31B)[iy_tay] , getTX(&Grid,G31B)[iy2_tay] ) )
			   -y* ( lerp( ay , getTX(&Grid,G22A)[iy_tay] , getTX(&Grid,G22A)[iy2_tay] ) ) ) ;
  const double dg1dy    = (lerp( ay , getTX(&Grid,G3Ady)[iy_tay] , getTX(&Grid,G3Ady)[iy2_tay] )
			   -lerp( ay , getTX(&Grid,G3Bdy)[iy_tay] , getTX(&Grid,G3Bdy)[iy2_tay] ) ) ;
  const double phi1     = ( lerp( ay , getTX(&Grid,G22A)[iy_tay] , getTX(&Grid,G22A)[iy2_tay] )
			    - lerp( ay , getTX(&Grid,G22B)[iy_tay] , getTX(&Grid,G22B)[iy2_tay] )/5. ) ;
  const double phi2     = ( lerp( ay , getTX(&Grid,G22B)[iy_tay] , getTX(&Grid,G22B)[iy2_tay] ) )*1.2 ;
  
  // various precomputations to make the loops faster
  const double ysq_4 = ysq/4. ;
  const double E = (ell3a + ell2) , D2 = (dell2adx + dell2ady) ,
    D1 = (dg1dy+dg1dx) , D12 = (dg1dx+dg2dx) , D3 = (y*ddg2dxdy-dg2dx) ;
  const double yvE[4]       = { yv[0]*E , yv[1]*E , yv[2]*E , yv[3]*E } ;
  const double yvell3a[4]   = { yv[0]*ell3a , yv[1]*ell3a , yv[2]*ell3a , yv[3]*ell3a } ;
  const double yhatD1[4]    = { yhat[0]*D1 , yhat[1]*D1 , yhat[2]*D1 , yhat[3]*D1 } ;
  const double yhatD2[4]    = { yhat[0]*D2 , yhat[1]*D2 , yhat[2]*D2 , yhat[3]*D2 } ;
  const double yhatD3[4]    = { yhat[0]*D3 , yhat[1]*D3 , yhat[2]*D3 , yhat[3]*D3 } ;
  const double yhatdg2dx[4] = { yhat[0]*dg2dx , yhat[1]*dg2dx , yhat[2]*dg2dx , yhat[3]*dg2dx } ;
  
  // point out the structures
  double *PtIII = (double*)tIII ;
  double *PtII  = (double*)tII ;
  double *PtI   = (double*)tI ;
  
  int alf , dta , bet ;
  for(alf=0;alf<4;alf++) {
    const double yalf2dx = yhat[alf]*(dell2adx) ;
    const double tIIsum  = dg0dx[alf]-0.5*yvell3a[alf]-ysq_4*yalf2dx;
    const double yD12    = yhat[alf]*D12 ;
    for(bet=0;bet<4;bet++) {
      const double dgSUM   = (dg0dx[bet]+dg0dy[bet]) ; 
      const double tIsum   = 2*(dlta[alf][bet]*phi1+yhat[alf]*yhat[bet]*phi2) ;
      const double tIIIsum = dgSUM-0.5*yvE[bet]-ysq_4*yhatD2[bet] ;
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
__device__ KQED_PRIVATE
int
Tabd_yeq0( const double xv[4] ,
	   const struct Grid_coeffs Grid ,
	   double tI[4][4][4] ,
	   double tII[4][4][4] ,
	   double tIII[4][4][4] )
{  
  const double xsq = SCALPROD(xv,xv);
  const double x = sqrt( fabs(xsq) ); 

  const int ix_tay = find_ind( Grid.XX , x , 0 , Grid.nstpx ) ;
  
  if( x > Grid.XX[ Grid.nstpx - 1 ] ) {
    // need to return here otherwise we will segfault
    return 1 ;
  }
  // we can set this now as we aren't at the upper edge
  const int ix2_tay = ix_tay+1 ;
  
  const double ax = ( Grid.XX[ix2_tay]-x)/( Grid.XX[ix2_tay]-Grid.XX[ix_tay]);
  const double xa = Grid.XX[ix_tay];
  const double xb = Grid.XX[ix2_tay];
  const double fx = 0.5*lerp( ax , getTY(&Grid,alpha0dx_0p)[ ix_tay ], getTY(&Grid,alpha0dx_0p)[ix2_tay] ) ;
  const double fy = ( lerp( ax , getTY(&Grid,alpha0_1p)[ix_tay] , getTY(&Grid,alpha0_1p)[ix2_tay] ) );

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

  const double ell1       = 4*( lerp( ax , getTY(&Grid,alpha1_0p)[ix_tay]/(xa4) ,
				      getTY(&Grid,alpha1_0p)[ix2_tay]/(xb4) ) )/3. ;
  const double ell3       = lerp( ax , getTY(&Grid,beta4_1p)[ix_tay]/(xa2) ,
				  getTY(&Grid,beta4_1p)[ix2_tay]/(xb2) ) ;
  
  const double dell1dx    = 4*lerp( ax , (getTY(&Grid,alpha1dx_0p)[ix_tay]
					  -4.*getTY(&Grid,alpha1_0p)[ix_tay]/xa)/xa4 ,
				    (getTY(&Grid,alpha1dx_0p)[ix2_tay]
				     -4.*getTY(&Grid,alpha1_0p)[ix2_tay]/xb)/xb4 )/3. ;
  
  const double dell1dycb  = 2*lerp( ax , (4.*getTY(&Grid,alpha1_1p)[ix_tay]/(3*xa)
					  - getTY(&Grid,beta4_1p)[ix_tay])/xa3 ,
				    (4.*getTY(&Grid,alpha1_1p)[ix2_tay]/(3*xb)
				     - getTY(&Grid,beta4_1p)[ix2_tay])/xb3 ) ;
  
  const double dg1dx      = lerp( ax , getTY(&Grid,alpha3dx_0p)[ix_tay] , getTY(&Grid,alpha3dx_0p)[ix2_tay] ) ;
  const double ddg1dxdx   = lerp( ax , getTY(&Grid,alpha3dxdx_0p)[ix_tay] , getTY(&Grid,alpha3dxdx_0p)[ix2_tay] ) ;
  const double dg2dx      = 2*lerp( ax , getTY(&Grid,beta2dx_1p)[ix_tay] , getTY(&Grid,beta2dx_1p)[ix2_tay] ) ;
  
  const double dg1dycb    = 2*lerp( ax , getTY(&Grid,alpha3_1p)[ix_tay]
				    - getTY(&Grid,beta2_1p)[ix_tay]/xa ,
				    getTY(&Grid,alpha3_1p)[ix2_tay]
				    - getTY(&Grid,beta2_1p)[ix2_tay]/xb ) ;
  
  const double ddg1dxdycb = 2*lerp( ax , getTY(&Grid,beta2_1p)[ix_tay]/(xa2)
				    - getTY(&Grid,beta2dx_1p)[ix_tay]/xa
				    + getTY(&Grid,alpha3dx_1p)[ix_tay] ,
				    getTY(&Grid,beta2_1p)[ix2_tay]/(xb2)
				    - getTY(&Grid,beta2dx_1p)[ix2_tay]/xb
				    + getTY(&Grid,alpha3dx_1p)[ix2_tay] ) ;

  // precomputations
  const double xsq_4 = xsq/4. , E = ell1+ell3 , D1 = dell1dx + dell1dycb
    , D2 = (ddg1dxdx+ddg1dxdycb) , D3 = (dg1dycb-dg2dx) , D4 = (dg1dx+dg1dycb) ;
  const double xvE[4]    = { xv[0]*E , xv[1]*E , xv[2]*E , xv[3]*E } ;
  const double xvell1[4] = { xv[0]*ell1 , xv[1]*ell1 , xv[2]*ell1 , xv[3]*ell1 } ;

  // point out the structures
  double *PtIII = (double*)tIII ;
  double *PtII  = (double*)tII ;
  double *PtI   = (double*)tI ;

  int alf , dta , bet ;
  for(alf=0;alf<4;alf++) {
    const double xhd1dx = xhat[alf]*(dell1dx) ;
    const double tIIsum = (dg0dx[alf]-0.5*xvell1[alf] ); 
    for(bet=0;bet<4;bet++) {
      const double Dg0SUM = (dg0dx[bet]+dg0dy[bet]-0.5*xvE[bet]) ;
      for(dta=0;dta<4;dta++) {
	const double xad = (xv[alf]*xv[dta]-xsq_4*dlta[alf][dta]);
	*PtIII = ((xad*xhat[bet]*D1 + dlta[alf][bet]*xvE[dta]+dlta[bet][dta]*xvE[alf]+Dg0SUM*dlta[alf][dta]) ); PtIII++ ;
	*PtII  = ((xv[bet]*xv[dta]-xsq_4*dlta[bet][dta])*xhd1dx+dlta[bet][alf]*xvell1[dta]+dlta[alf][dta]*xvell1[bet]+dlta[bet][dta]*(tIIsum)); PtII++ ;
	*PtI   = ((dlta[alf][dta]*xhat[bet]+dlta[bet][dta]*xhat[alf]+dlta[alf][bet]*xhat[dta]-xhat[alf]*xhat[bet]*xhat[dta])*D4+xhat[alf]*( xhat[bet]*xv[dta]*D2 - dlta[bet][dta]*D3 ));PtI++ ;	
      }
    }
  }
  
  return 0 ;
}
