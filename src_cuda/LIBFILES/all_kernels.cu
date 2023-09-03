/**
   @file all_kernels.c
   @brief compute all kernels L^{0/1/2/3} for (x,y) and (y,x) arguements

   Also contains the average function for all kernels
 */
#include "KQED.h"

#include "kernels.h"

// from linear idx in 0->384 return a value where mu and lambda indices are swapped
__device__
static inline size_t
i_to_mulam( const size_t idx )
{
  const size_t l[4] = { idx/64 , (idx/16)&3 , (idx/4)&3 , idx&3 } ;
  return l[1] + 4*(l[2]+4*(l[3]+4*l[0]) ) ;
}

// computes *K = 0.5*( *K + *W ) for all kernels in K and W
__device__
void
average_all_kernels( struct Kernels *K ,
		     const struct Kernels W )
{
  // point out elements of K
  double *KL0 = (double*)K->L0 , *KL1 = (double*)K->L1 ;
  double *KL2 = (double*)K->L2 , *KL3 = (double*)K->L3 ;
  
  // const point out elements of W
  const double *WL0 = (const double*)W.L0 ;
  const double *WL1 = (const double*)W.L1 ;
  const double *WL2 = (const double*)W.L2 ;
  const double *WL3 = (const double*)W.L3 ;

  // do the average
  size_t i ;
  for( i = 0 ; i < 384 ; i++ ) {
    *KL0 = 0.5*( *KL0 + *WL0 ) ; KL0++ ; WL0++ ;
    *KL1 = 0.5*( *KL1 + *WL1 ) ; KL1++ ; WL1++ ;
    *KL2 = 0.5*( *KL2 + *WL2 ) ; KL2++ ; WL2++ ;
    *KL3 = 0.5*( *KL3 + *WL3 ) ; KL3++ ; WL3++ ;
  }
  return ;
}

// computes *K = 0.5*( *K + *W ) for all kernels in K and W
__device__
void
average_all_QED_kernels( struct QED_Kernels *K ,
			 const struct QED_Kernels W )
{
  // point out elements of K
  double *KL0xy = (double*)K->L0.xy ;
  double *KL0yx = (double*)K->L0.yx ;
  double *KL1xy = (double*)K->L1.xy ;
  double *KL1yx = (double*)K->L1.yx ;
  double *KL2xy = (double*)K->L2.xy ;
  double *KL2yx = (double*)K->L2.yx ;
  double *KL3xy = (double*)K->L3.xy ;
  double *KL3yx = (double*)K->L3.yx ;
  
  // const point out elements of W
  const double *WL0xy = (const double*)W.L0.xy ;
  const double *WL0yx = (const double*)W.L0.yx ;
  const double *WL1xy = (const double*)W.L1.xy ;
  const double *WL1yx = (const double*)W.L1.yx ;
  const double *WL2xy = (const double*)W.L2.xy ;
  const double *WL2yx = (const double*)W.L2.yx ;
  const double *WL3xy = (const double*)W.L3.xy ;
  const double *WL3yx = (const double*)W.L3.yx ;

  // do the average
  size_t i ;
  for( i = 0 ; i < 384 ; i++ ) {
    *KL0xy = 0.5*( *KL0xy + *WL0xy ) ; KL0xy++ ; WL0xy++ ;
    *KL0yx = 0.5*( *KL0yx + *WL0yx ) ; KL0yx++ ; WL0yx++ ;
    *KL1xy = 0.5*( *KL1xy + *WL1xy ) ; KL1xy++ ; WL1xy++ ;
    *KL1yx = 0.5*( *KL1yx + *WL1yx ) ; KL1yx++ ; WL1yx++ ;
    *KL2xy = 0.5*( *KL2xy + *WL2xy ) ; KL2xy++ ; WL2xy++ ;
    *KL2yx = 0.5*( *KL2yx + *WL2yx ) ; KL2yx++ ; WL2yx++ ;
    *KL3xy = 0.5*( *KL3xy + *WL3xy ) ; KL3xy++ ; WL3xy++ ;
    *KL3yx = 0.5*( *KL3yx + *WL3yx ) ; KL3yx++ ; WL3yx++ ;
  }
  return ;
}

// computes L_{rhosig;mu,nu,lambda}(x,y) and L_{rhosig;mu,nu,lambda}(y,x)
__device__
void
compute_all_kernels( const double xv[4] ,
		     const double yv[4] ,
		     const struct QED_kernel_temps t ,
		     struct QED_Kernels *K )
{
  // we only need 6 kernels to compute everything
  // K[0].Lxy = L(x,y) K[0].Lyx = L(y,x)
  QED_kernel_L0( xv , yv , t , K -> L0.xy ) ;
  QED_kernel_L0( yv , xv , t , K -> L0.yx ) ;

  // and these temporaries: L(x,0) L(y,0) L(0,x) L(0,y)
  double Lx0[6][4][4][4] , L0x[6][4][4][4] ;
  double Ly0[6][4][4][4] , L0y[6][4][4][4] ;
  const double zero[4] = { 0. , 0. , 0. , 0. } ;

  QED_kernel_L0( xv , zero , t , Lx0 ) ;
  QED_kernel_L0( zero , xv , t , L0x ) ;

  QED_kernel_L0( yv , zero , t , Ly0 ) ;
  QED_kernel_L0( zero , yv , t , L0y ) ;

  // point out the data structures we want
  const double *L0xy   = (const double*)K->L0.xy ;
  const double *L0yx   = (const double*)K->L0.yx ;
  const double *pLx0   = (const double*)Lx0 ;
  const double *pLy0   = (const double*)Ly0 ;
  const double *pLx0ml = (const double*)Lx0 ;
  const double *pLy0ml = (const double*)Ly0 ;
  const double *pL0x   = (const double*)L0x ;
  const double *pL0y   = (const double*)L0y ;

  // write into these structures
  double *L1xy = (double*)K->L1.xy ;
  double *L2xy = (double*)K->L2.xy ;
  double *L3xy = (double*)K->L3.xy ;
  double *L1yx = (double*)K->L1.yx ;
  double *L2yx = (double*)K->L2.yx ;
  double *L3yx = (double*)K->L3.yx ;
  
  size_t i ;
  for( i = 0 ; i < 384 ; i++ ) {

    const size_t idx = i_to_mulam(i) ;
    const double Lx0mulam = *(pLx0ml + idx) ;
    const double Ly0mulam = *(pLy0ml + idx) ;
    const double f1 = ( Lx0mulam + Ly0mulam )/2. ;
    const double f2 = ( *pLx0 + *pL0y ) ;
    const double f4 = ( *pL0x + *pLy0 ) ;
    const double f3 = ( *pL0x - *pL0y ) ;

    // L^(1)(x,y)
    *L1xy = *L0xy + f1 ;
    // L^(2)(x,y)
    *L2xy = *L0xy - f2 ;
    // L^3(x,y)
    *L3xy = *L0xy + Lx0mulam + f3 ;
    // L^1(y,x)
    *L1yx = *L0yx + f1 ;
    // L^2(y,x)
    *L2yx = *L0yx - f4 ;
    // L^3(y,x)
    *L3yx = *L0yx + Ly0mulam - f3 ;
    
    // increment pointers
    L0xy++ ; L1xy++ ; L2xy++ ; L3xy++ ;
    L0yx++ ; L1yx++ ; L2yx++ ; L3yx++ ;
    pLx0 ++ ; pLy0++ ; pL0x++ ; pL0y++ ;
  }
  
  return ;
}

// copy for the XY
#define inlineXY(a)							\
  ( + eX[a]*(L0y[rhosig][mu][nu][lambda] -				\
	     Mx[mu][a]*(xv[0]*L0y[rhosig][0][nu][lambda] +		\
			xv[1]*L0y[rhosig][1][nu][lambda] +		\
			xv[2]*L0y[rhosig][2][nu][lambda] +		\
			xv[3]*L0y[rhosig][3][nu][lambda] ) )		\
    + eY[a]*(Lx0[rhosig][mu][nu][lambda] -				\
	     My[nu][a]*(yv[0]*Lx0[rhosig][mu][0][lambda] +		\
			yv[1]*Lx0[rhosig][mu][1][lambda] +		\
			yv[2]*Lx0[rhosig][mu][2][lambda] +		\
			yv[3]*Lx0[rhosig][mu][3][lambda] ) ) )

// copy for the XY
#define inlineXMY(a)							\
  ( - eX[a]*(L0y[rhosig][mu][nu][lambda] -				\
	     Mx[mu][a]*(xv[0]*L0y[rhosig][0][nu][lambda] +		\
			xv[1]*L0y[rhosig][1][nu][lambda] +		\
			xv[2]*L0y[rhosig][2][nu][lambda] +		\
			xv[3]*L0y[rhosig][3][nu][lambda] ) )		\
    + eY[a]*(Lx0[rhosig][mu][nu][lambda] -				\
	     My[nu][a]*(yv[0]*Lx0[rhosig][mu][0][lambda] +		\
			yv[1]*Lx0[rhosig][mu][1][lambda] +		\
			yv[2]*Lx0[rhosig][mu][2][lambda] +		\
			yv[3]*Lx0[rhosig][mu][3][lambda] ) ) )

// copy for the YX
#define inlineYX(a)				\
  ( + eY[a]*(L0x[rhosig][mu][nu][lambda] -				\
	   My[mu][a]*(yv[0]*L0x[rhosig][0][nu][lambda] +		\
		      yv[1]*L0x[rhosig][1][nu][lambda] +		\
		      yv[2]*L0x[rhosig][2][nu][lambda] +		\
		      yv[3]*L0x[rhosig][3][nu][lambda] ) )		\
    + eX[a]*(Ly0[rhosig][mu][nu][lambda] -				\
	     Mx[nu][a]*(xv[0]*Ly0[rhosig][mu][0][lambda] +		\
			xv[1]*Ly0[rhosig][mu][1][lambda] +		\
			xv[2]*Ly0[rhosig][mu][2][lambda] +		\
			xv[3]*Ly0[rhosig][mu][3][lambda] ) ) )

// copy for the YX
#define inlineYX2(a)							\
  ( + eY[a]*(L0x[rhosig][nu][mu][lambda] -				\
	   My[nu][a]*(yv[0]*L0x[rhosig][0][mu][lambda] +		\
		      yv[1]*L0x[rhosig][1][mu][lambda] +		\
		      yv[2]*L0x[rhosig][2][mu][lambda] +		\
		      yv[3]*L0x[rhosig][3][mu][lambda] ) )		\
    + eX[a]*(Ly0[rhosig][nu][mu][lambda] -				\
	     Mx[mu][a]*(xv[0]*Ly0[rhosig][nu][0][lambda] +		\
			xv[1]*Ly0[rhosig][nu][1][lambda] +		\
			xv[2]*Ly0[rhosig][nu][2][lambda] +		\
			xv[3]*Ly0[rhosig][nu][3][lambda] ) ) )

// computes L_{rhosig;mu,nu,lambda}(x,y) and L_{rhosig;mu,nu,lambda}(y,x)
__device__
void
compute_all_Mkernels( const double M[4] ,
		      const double xv[4] ,
		      const double yv[4] ,
		      const struct QED_kernel_temps t ,
		      struct QED_Kernels *K )
{
  // we only need 6 kernels to compute everything
  double Lxy[6][4][4][4] , Lyx[6][4][4][4] ;
  QED_kernel_L0( xv , yv , t , Lxy ) ;
  QED_kernel_L0( yv , xv , t , Lyx ) ;

  // and these temporaries: L(x,0) L(y,0) L(0,x) L(0,y)
  double Lx0[6][4][4][4] , L0x[6][4][4][4] ;
  double Ly0[6][4][4][4] , L0y[6][4][4][4] ;
  const double zero[4] = { 0. , 0. , 0. , 0. } ;

  QED_kernel_L0( xv , zero , t , Lx0 ) ;
  QED_kernel_L0( zero , xv , t , L0x ) ;

  QED_kernel_L0( yv , zero , t , Ly0 ) ;
  QED_kernel_L0( zero , yv , t , L0y ) ;

    // precompute x^2, y^2
  const double xsq = xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]+xv[3]*xv[3] ;
  const double ysq = yv[0]*yv[0]+yv[1]*yv[1]+yv[2]*yv[2]+yv[3]*yv[3] ;

  // precompute the gaussians
  const double eX[4] = { exp(-M[0]*xsq/2.) , exp(-M[1]*xsq/2.) ,
			 exp(-M[2]*xsq/2.) , exp(-M[3]*xsq/2.) } ;
  const double eY[4] = { exp(-M[0]*ysq/2.) , exp(-M[1]*ysq/2.) ,
			 exp(-M[2]*ysq/2.) , exp(-M[3]*ysq/2.) } ;

  // look-up table precomputations
  const double Mx[4][4] = { { M[0]*xv[0] , M[1]*xv[0] , M[2]*xv[0] , M[3]*xv[0] } ,
			    { M[0]*xv[1] , M[1]*xv[1] , M[2]*xv[1] , M[3]*xv[1] } ,
			    { M[0]*xv[2] , M[1]*xv[2] , M[2]*xv[2] , M[3]*xv[2] } ,
			    { M[0]*xv[3] , M[1]*xv[3] , M[2]*xv[3] , M[3]*xv[3] } } ;
  const double My[4][4] = { { M[0]*yv[0] , M[1]*yv[0] , M[2]*yv[0] , M[3]*yv[0] } ,
			    { M[0]*yv[1] , M[1]*yv[1] , M[2]*yv[1] , M[3]*yv[1] } ,
			    { M[0]*yv[2] , M[1]*yv[2] , M[2]*yv[2] , M[3]*yv[2] } ,
			    { M[0]*yv[3] , M[1]*yv[3] , M[2]*yv[3] , M[3]*yv[3] } } ;

  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // M0
	  K->L0.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(0) ;
	  K->L0.yx[rhosig][mu][nu][lambda] =
	    Lyx[rhosig][mu][nu][lambda] - inlineYX(0) ;
	  // M1
	  K->L1.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(1) ;
	  K->L1.yx[rhosig][mu][nu][lambda] =
	    Lyx[rhosig][mu][nu][lambda] - inlineYX(1) ;
	  // M2
	  K->L2.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(2) ;
	  K->L2.yx[rhosig][mu][nu][lambda] =
	    Lyx[rhosig][mu][nu][lambda] - inlineYX(2) ;
	  // M3
	  K->L3.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(3) ;
	  K->L3.yx[rhosig][mu][nu][lambda] =
	    Lyx[rhosig][mu][nu][lambda] - inlineYX(3) ;
	}
      }
    }
  }
  return ;
}

// computes L_{rhosig;mu,nu,lambda}(x,y) + L_{rhosig;nu,mu,lambda}(y,x) in Lxy
// and puts L_{rhosig;mu,nu,lambda}(x,-y) in Lyx
__device__
void
compute_all_Mkernels_v2( const double M[4] ,
			 const double xv[4] ,
			 const double yv[4] ,
			 const struct QED_kernel_temps t ,
			 struct QED_Kernels *K )
{
  // we only need 6 kernels to compute everything
  double Lxy[6][4][4][4] , Lyx[6][4][4][4] , Lxmy[6][4][4][4] ;

  const double myv[4] = { -yv[0] , -yv[1] , -yv[2] , -yv[3] } ;
  QED_kernel_L0( xv , yv  , t , Lxy  ) ;
  QED_kernel_L0( xv , myv , t , Lxmy ) ;
  QED_kernel_L0( yv , xv  , t , Lyx  ) ;

  // and these temporaries: L(x,0) L(y,0) L(0,x) L(0,y)
  double Lx0[6][4][4][4] , L0x[6][4][4][4] ;
  double Ly0[6][4][4][4] , L0y[6][4][4][4] ;
  const double zero[4] = { 0. , 0. , 0. , 0. } ;

  QED_kernel_L0( xv , zero , t , Lx0 ) ;
  QED_kernel_L0( zero , xv , t , L0x ) ;

  QED_kernel_L0( yv , zero , t , Ly0 ) ;
  QED_kernel_L0( zero , yv , t , L0y ) ;

    // precompute x^2, y^2
  const double xsq = xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]+xv[3]*xv[3] ;
  const double ysq = yv[0]*yv[0]+yv[1]*yv[1]+yv[2]*yv[2]+yv[3]*yv[3] ;

  // precompute the gaussians
  const double eX[4] = { exp(-M[0]*xsq/2.) , exp(-M[1]*xsq/2.) ,
			 exp(-M[2]*xsq/2.) , exp(-M[3]*xsq/2.) } ;
  const double eY[4] = { exp(-M[0]*ysq/2.) , exp(-M[1]*ysq/2.) ,
			 exp(-M[2]*ysq/2.) , exp(-M[3]*ysq/2.) } ;

  // look-up table precomputations
  const double Mx[4][4] = { { M[0]*xv[0] , M[1]*xv[0] , M[2]*xv[0] , M[3]*xv[0] } ,
			    { M[0]*xv[1] , M[1]*xv[1] , M[2]*xv[1] , M[3]*xv[1] } ,
			    { M[0]*xv[2] , M[1]*xv[2] , M[2]*xv[2] , M[3]*xv[2] } ,
			    { M[0]*xv[3] , M[1]*xv[3] , M[2]*xv[3] , M[3]*xv[3] } } ;
  const double My[4][4] = { { M[0]*yv[0] , M[1]*yv[0] , M[2]*yv[0] , M[3]*yv[0] } ,
			    { M[0]*yv[1] , M[1]*yv[1] , M[2]*yv[1] , M[3]*yv[1] } ,
			    { M[0]*yv[2] , M[1]*yv[2] , M[2]*yv[2] , M[3]*yv[2] } ,
			    { M[0]*yv[3] , M[1]*yv[3] , M[2]*yv[3] , M[3]*yv[3] } } ;

  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // M0
	  K->L0.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(0) +
	    Lyx[rhosig][nu][mu][lambda] - inlineYX2(0) ;
	  // M1
	  K->L1.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(1) +
	    Lyx[rhosig][nu][mu][lambda] - inlineYX2(1) ;
	  // M2
	  K->L2.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(2) +
	    Lyx[rhosig][nu][mu][lambda] - inlineYX2(2) ;
	  // M3
	  K->L3.xy[rhosig][mu][nu][lambda] =
	    Lxy[rhosig][mu][nu][lambda] - inlineXY(3) +
	    Lyx[rhosig][nu][mu][lambda] - inlineYX2(3) ;
	  // and the -y terms
	  K->L0.yx[rhosig][mu][nu][lambda] =
	    Lxmy[rhosig][mu][nu][lambda] - inlineXMY(0) ;
	  K->L1.yx[rhosig][mu][nu][lambda] =
	    Lxmy[rhosig][mu][nu][lambda] - inlineXMY(1) ;
	  K->L2.yx[rhosig][mu][nu][lambda] =
	    Lxmy[rhosig][mu][nu][lambda] - inlineXMY(2) ;
	  K->L3.yx[rhosig][mu][nu][lambda] =
	    Lxmy[rhosig][mu][nu][lambda] - inlineXMY(3) ;
	}
      }
    }
  }
  return ;
}

// swaps L_{..mu,nu..}(y,x) and  L_{..nu,mu..}(y,x) for all kernels
__device__
void
swap_munu_Lyx( struct QED_Kernels *K )
{
  size_t rhosig , mu , nu , lambda ;
  double tmp ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = mu+1 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // L0 kernel
	  tmp = K -> L0.yx[rhosig][nu][mu][lambda] ;
	  K -> L0.yx[rhosig][nu][mu][lambda] =
	    K -> L0.yx[rhosig][mu][nu][lambda] ;
	  K -> L0.yx[rhosig][mu][nu][lambda] = tmp ;
	  // L1 kernel
	  tmp = K -> L1.yx[rhosig][nu][mu][lambda] ;
	  K -> L1.yx[rhosig][nu][mu][lambda] =
	    K -> L1.yx[rhosig][mu][nu][lambda] ;
	  K -> L1.yx[rhosig][mu][nu][lambda] = tmp ;
	  // L2 kernel
	  tmp = K -> L2.yx[rhosig][nu][mu][lambda] ;
	  K -> L2.yx[rhosig][nu][mu][lambda] =
	    K -> L2.yx[rhosig][mu][nu][lambda] ;
	  K -> L2.yx[rhosig][mu][nu][lambda] = tmp ;
	  // L3 kernel
	  tmp = K -> L3.yx[rhosig][nu][mu][lambda] ;
	  K -> L3.yx[rhosig][nu][mu][lambda] =
	    K -> L3.yx[rhosig][mu][nu][lambda] ;
	  K -> L3.yx[rhosig][mu][nu][lambda] = tmp ;
	}
      }
    }
  }	
  
  return ;
}
