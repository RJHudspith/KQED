/**
   @file con_kernel.c
   @brief computes the kernel needed for the connected contribution

   ( dropping rho,sigma arguments)
   
   A = L_{\mu\nu\lambda}(x,y)
   B = L_{\lambda\nu\mu}(x,x-y)
   C = L_{\mu\nu\lambda}(x,y) + L_{\nu\mu\lambda}(y,x) - L_{\lambda\nu\mu}(x,x-y)
*/
#include "KQED.h"

#include "kernels.h"

// computes *K = 0.5*( *K + *W ) for all con kernels in K and W
__device__
void
average_con_kernels( struct Kernels *K ,
		     const struct Kernels *W )
{
  // point out elements of K
  double *KL0A = (double*)K[0].L0 , *KL0B = (double*)K[1].L0 ;
  double *KL0C = (double*)K[2].L0 , *KL1A = (double*)K[0].L1 ;
  double *KL1B = (double*)K[1].L1 , *KL1C = (double*)K[2].L1 ;
  double *KL2A = (double*)K[0].L2 , *KL2B = (double*)K[1].L2 ;
  double *KL2C = (double*)K[2].L2 , *KL3A = (double*)K[0].L3 ;
  double *KL3B = (double*)K[1].L3 , *KL3C = (double*)K[2].L3 ;
  
  // const point out elements of W
  const double *WL0A = (const double*)W[0].L0 ;
  const double *WL1A = (const double*)W[0].L1 ;
  const double *WL2A = (const double*)W[0].L2 ;
  const double *WL3A = (const double*)W[0].L3 ;
  const double *WL3B = (const double*)W[1].L3 ;
  const double *WL0B = (const double*)W[1].L0 ;
  const double *WL1B = (const double*)W[1].L1 ;
  const double *WL2B = (const double*)W[1].L2 ; 
  const double *WL0C = (const double*)W[2].L0 ;
  const double *WL1C = (const double*)W[2].L1 ;
  const double *WL2C = (const double*)W[2].L2 ;
  const double *WL3C = (const double*)W[2].L3 ;
  
  // do the average
  size_t i ;
  for( i = 0 ; i < 384 ; i++ ) {
    *KL0A = 0.5*( *KL0A + *WL0A ) ; KL0A++ ; WL0A++ ;
    *KL0B = 0.5*( *KL0B + *WL0B ) ; KL0B++ ; WL0B++ ;
    *KL0C = 0.5*( *KL0C + *WL0C ) ; KL0C++ ; WL0C++ ;

    *KL1A = 0.5*( *KL1A + *WL1A ) ; KL1A++ ; WL1A++ ;
    *KL1B = 0.5*( *KL1B + *WL1B ) ; KL1B++ ; WL1B++ ;
    *KL1C = 0.5*( *KL1C + *WL1C ) ; KL1C++ ; WL1C++ ;

    *KL2A = 0.5*( *KL2A + *WL2A ) ; KL2A++ ; WL2A++ ;
    *KL2B = 0.5*( *KL2B + *WL2B ) ; KL2B++ ; WL2B++ ;
    *KL2C = 0.5*( *KL2C + *WL2C ) ; KL2C++ ; WL2C++ ;
    
    *KL3A = 0.5*( *KL3A + *WL3A ) ; KL3A++ ; WL3A++ ;
    *KL3B = 0.5*( *KL3B + *WL3B ) ; KL3B++ ; WL3B++ ;
    *KL3C = 0.5*( *KL3C + *WL3C ) ; KL3C++ ; WL3C++ ;
  }
  return ;
}	 

__device__
void
compute_con_kernels_v2( const double xv[4] ,
			const double yv[4] ,
			const struct QED_kernel_temps t ,
			struct Kernels *K )
{
  const double xmyv[4] = { xv[0]-yv[0] , xv[1]-yv[1] ,
			   xv[2]-yv[2] , xv[3]-yv[3] } ;
  const double zero[4] = { 0. , 0. , 0. , 0. } ;
  double Kx_xmy[6][4][4][4] , Ky_x[6][4][4][4] ;
  QED_kernel_L0( xv , yv , t , K[0].L0 ) ;
  QED_kernel_L0( yv , xv , t , Ky_x ) ;
  QED_kernel_L0( xv , xmyv , t , Kx_xmy ) ;

  // and then the ones with a zero
  double Kx_0[6][4][4][4]   , K0_x[6][4][4][4] ;
  double Ky_0[6][4][4][4]   , K0_y[6][4][4][4] ;
  double Kxmy_0[6][4][4][4] , K0_xmy[6][4][4][4] ;
  QED_kernel_L0( xv , zero , t , Kx_0 ) ;
  QED_kernel_L0( zero , xv , t , K0_x ) ;
  QED_kernel_L0( yv , zero , t , Ky_0 ) ;
  QED_kernel_L0( zero , yv , t , K0_y ) ;
  QED_kernel_L0( xmyv , zero , t , Kxmy_0 ) ;
  QED_kernel_L0( zero , xmyv , t , K0_xmy ) ;

  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // L0 kernel
	  // x,x-y
	  K[1].L0[rhosig][mu][nu][lambda] =
	    Kx_xmy[rhosig][lambda][nu][mu] ; 
	  K[2].L0[rhosig][mu][nu][lambda] =
	    K[0].L0[rhosig][mu][nu][lambda]
	      +Ky_x[rhosig][nu][mu][lambda]
	    -Kx_xmy[rhosig][lambda][nu][mu] ;
	  // L1 kernel
	  K[0].L1[rhosig][mu][nu][lambda] = 
	    K[0].L0[rhosig][mu][nu][lambda] +
	    ( +Kx_0[rhosig][lambda][nu][mu]
	      +Ky_0[rhosig][lambda][nu][mu] )/2 ;
	  K[1].L1[rhosig][mu][nu][lambda] = 
	      K[1].L0[rhosig][mu][nu][lambda] +
	    (   +Kx_0[rhosig][mu][nu][lambda]
	      +Kxmy_0[rhosig][mu][nu][lambda] )/2 ;
	  K[2].L1[rhosig][mu][nu][lambda] = 
	    K[0].L1[rhosig][mu][nu][lambda] -
	    K[1].L1[rhosig][mu][nu][lambda] +
	       Ky_x[rhosig][nu][mu][lambda] +
	    ( +Ky_0[rhosig][lambda][mu][nu]
	      +Kx_0[rhosig][lambda][mu][nu] )/2 ;
	  // L2 kernel
	  K[0].L2[rhosig][mu][nu][lambda] = 
	    K[0].L0[rhosig][mu][nu][lambda] -
	    ( +K0_y[rhosig][mu][nu][lambda]
	      +Kx_0[rhosig][mu][nu][lambda] ) ;
	  K[1].L2[rhosig][mu][nu][lambda] = 
	      K[1].L0[rhosig][mu][nu][lambda] -
	    ( +K0_xmy[rhosig][lambda][nu][mu]
	        +Kx_0[rhosig][lambda][nu][mu] ) ;
	  K[2].L2[rhosig][mu][nu][lambda] =
	    K[0].L2[rhosig][mu][nu][lambda] -
	    K[1].L2[rhosig][mu][nu][lambda] +
	       Ky_x[rhosig][nu][mu][lambda] -
	    ( +K0_x[rhosig][nu][mu][lambda]
	      +Ky_0[rhosig][nu][mu][lambda] ) ;
	  // L3 kernel
	  K[0].L3[rhosig][mu][nu][lambda] = 
	    K[0].L0[rhosig][mu][nu][lambda] +
	    ( +Kx_0[rhosig][lambda][nu][mu]
	      +K0_x[rhosig][mu][nu][lambda]
	      -K0_y[rhosig][mu][nu][lambda] ) ;
	  K[1].L3[rhosig][mu][nu][lambda] = 
	      K[1].L0[rhosig][mu][nu][lambda] +
	    (   +Kx_0[rhosig][mu][nu][lambda]
	        +K0_x[rhosig][lambda][nu][mu]
	      -K0_xmy[rhosig][lambda][nu][mu] ) ; 
	  K[2].L3[rhosig][mu][nu][lambda] =
	    K[0].L3[rhosig][mu][nu][lambda] -
	    K[1].L3[rhosig][mu][nu][lambda] +
	       Ky_x[rhosig][nu][mu][lambda] +
	    ( +Ky_0[rhosig][lambda][mu][nu]
	      +K0_y[rhosig][nu][mu][lambda]
	      -K0_x[rhosig][nu][mu][lambda] ) ; 	  
	}
      }
    }
  }
  return ;
}

__device__
void
compute_con_kernels( const double xv[4] ,
		     const double yv[4] ,
		     const struct QED_kernel_temps t ,
		     struct QED_Kernels *K )
{
  const double xmyv[4] = { xv[0]-yv[0] , xv[1]-yv[1] ,
			   xv[2]-yv[2] , xv[3]-yv[3] } ;
  const double zero[4] = { 0. , 0. , 0. , 0. } ;
  double Kx_xmy[6][4][4][4] , Ky_x[6][4][4][4] ;
  QED_kernel_L0( xv , yv , t , K->L0.xy ) ;
  QED_kernel_L0( yv , xv , t , Ky_x ) ;
  QED_kernel_L0( xv , xmyv , t , Kx_xmy ) ;

  // and then the ones with a zero
  double Kx_0[6][4][4][4]   , K0_x[6][4][4][4] ;
  double Ky_0[6][4][4][4]   , K0_y[6][4][4][4] ;
  double Kxmy_0[6][4][4][4] , K0_xmy[6][4][4][4] ;
  QED_kernel_L0( xv , zero , t , Kx_0 ) ;
  QED_kernel_L0( zero , xv , t , K0_x ) ;
  QED_kernel_L0( yv , zero , t , Ky_0 ) ;
  QED_kernel_L0( zero , yv , t , K0_y ) ;
  QED_kernel_L0( xmyv , zero , t , Kxmy_0 ) ;
  QED_kernel_L0( zero , xmyv , t , K0_xmy ) ;

  size_t rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  // L0 kernel
	  K->L0.xy[rhosig][mu][nu][lambda] =
	    ( +K->L0.xy[rhosig][mu][nu][lambda]
	      +Ky_x[rhosig][nu][mu][lambda]
	      -Kx_xmy[rhosig][lambda][nu][mu] ) ;
	  K->L0.yx[rhosig][mu][nu][lambda] =
	    Kx_xmy[rhosig][lambda][nu][mu] ; 

	  // L1 kernel
	  K->L1.xy[rhosig][mu][nu][lambda] =
	    K->L0.xy[rhosig][mu][nu][lambda] +
	    ( +Kx_0[rhosig][lambda][nu][mu]
	      +Ky_0[rhosig][lambda][nu][mu]
	      +Kx_0[rhosig][lambda][mu][nu]
	      +Ky_0[rhosig][lambda][mu][nu]
	      -Kx_0[rhosig][mu][nu][lambda]
	      -Kxmy_0[rhosig][mu][nu][lambda] )/2 ;
	  K->L1.yx[rhosig][mu][nu][lambda] = 
	    K->L0.yx[rhosig][mu][nu][lambda] +
	    ( +Kx_0[rhosig][mu][nu][lambda]
	      +Kxmy_0[rhosig][mu][nu][lambda] )/2 ;
	  
	  // L2 kernel
	  K->L2.xy[rhosig][mu][nu][lambda] =
	    K->L0.xy[rhosig][mu][nu][lambda] -
	    ( +K0_y[rhosig][mu][nu][lambda]
	      +Kx_0[rhosig][mu][nu][lambda]
	      +K0_x[rhosig][nu][mu][lambda]
	      +Ky_0[rhosig][nu][mu][lambda]
	      -K0_xmy[rhosig][lambda][nu][mu]
	      -Kx_0[rhosig][lambda][nu][mu] ) ;
	  K->L2.yx[rhosig][mu][nu][lambda] = 
	    K->L0.yx[rhosig][mu][nu][lambda] -
	    ( +Kx_0[rhosig][lambda][nu][mu]
	      +K0_xmy[rhosig][lambda][nu][mu] ) ;

	  // L3 kernel
	  K->L3.xy[rhosig][mu][nu][lambda] =
	    K->L2.xy[rhosig][mu][nu][lambda] +
	    ( K0_x[rhosig][mu][nu][lambda] +
	      K0_y[rhosig][nu][mu][lambda] +
	      Ky_0[rhosig][lambda][mu][nu] -
	      K0_x[rhosig][lambda][nu][mu] +
	      Ky_0[rhosig][nu][mu][lambda] ) ;
	  K->L3.yx[rhosig][mu][nu][lambda] = 
	    K->L0.yx[rhosig][mu][nu][lambda] +
	    ( +Kx_0[rhosig][mu][nu][lambda]
	      +K0_x[rhosig][lambda][nu][mu]
	      -K0_xmy[rhosig][lambda][nu][mu] ) ;
	}
      }
    }
  }
	  
  return ;
}

// hmmm apparently this is pretty slow, perhaps I need to do something
// drastic about it
#define UNROLL_NM(Nm)							\
  LMxy[Nm] = Kxy[rhosig][mu][nu][lambda]				\
    -*eXP*( K0y[rhosig][mu][nu][lambda] - Mx[mu][Nm]*Pf1 )		\
    -*eYP*( Kx0[rhosig][mu][nu][lambda] - My[nu][Nm]*Pf2[lambda] )  ;	\
  LMyx[Nm] = Kyx[rhosig][nu][mu][lambda]				\
    -*eYP*( K0x[rhosig][nu][mu][lambda]-My[nu][Nm]*Pf3[lambda] )	\
    -*eXP*( Ky0[rhosig][nu][mu][lambda]-Mx[mu][Nm]*Pf4 )  ;		\
  LMxxmy[Nm] = Kxxmy[rhosig][lambda][nu][mu]				\
    -*eXP*( K0xmy[rhosig][lambda][nu][mu]-Mx[lambda][Nm]*Pf5 )		\
    -*eXMYP*( Kx0[rhosig][lambda][nu][mu]-Mxmy[nu][Nm]*Pf6[lambda] )  ; \
  eXP++ ; eYP++ ; eXMYP++ ;

// compute 4 L2 kernels with the M-factor included
__device__
void
compute_con_kernelsM_L2( const double M[4] ,
			 const double xv[4] ,
			 const double yv[4] ,
			 const struct QED_kernel_temps t ,
			 struct Kernels *K )
{
  const double xmyv[4] = { xv[0]-yv[0] , xv[1]-yv[1] , xv[2]-yv[2] , xv[3]-yv[3] } ;
  
  double Kxy[6][4][4][4] KQED_ALIGN , Kxxmy[6][4][4][4] KQED_ALIGN ;
  double Kyx[6][4][4][4] KQED_ALIGN ;
  QED_kernel_L0( xv , yv   , t , Kxy ) ;
  QED_kernel_L0( yv , xv   , t , Kyx ) ;
  QED_kernel_L0( xv , xmyv , t , Kxxmy ) ;

  double Kx0[6][4][4][4] KQED_ALIGN , Ky0[6][4][4][4] KQED_ALIGN ;
  double K0y[6][4][4][4] KQED_ALIGN , K0x[6][4][4][4] KQED_ALIGN ;
  double K0xmy[6][4][4][4] KQED_ALIGN ;
  const double zero[4] = { 0. , 0. , 0. , 0. } ;
  QED_kernel_L0( xv , zero   , t , Kx0   ) ;
  QED_kernel_L0( zero , xv   , t , K0x   ) ;
  QED_kernel_L0( yv , zero   , t , Ky0   ) ;
  QED_kernel_L0( zero , yv   , t , K0y   ) ;
  QED_kernel_L0( zero , xmyv , t , K0xmy ) ;

  // precompute x^2, y^2, (x-y)^2
  const double xsq = xv[0]*xv[0]+xv[1]*xv[1]+xv[2]*xv[2]+xv[3]*xv[3] ;
  const double ysq = yv[0]*yv[0]+yv[1]*yv[1]+yv[2]*yv[2]+yv[3]*yv[3] ;
  const double xmysq = xmyv[0]*xmyv[0]+xmyv[1]*xmyv[1]+xmyv[2]*xmyv[2]+xmyv[3]*xmyv[3] ;

  // precompute the gaussians
  const double eX[4] KQED_ALIGN = { exp(-M[0]*xsq/2.) , exp(-M[1]*xsq/2.) ,
				    exp(-M[2]*xsq/2.) , exp(-M[3]*xsq/2.) } ;
  const double eY[4] KQED_ALIGN = { exp(-M[0]*ysq/2.) , exp(-M[1]*ysq/2.) ,
				    exp(-M[2]*ysq/2.) , exp(-M[3]*ysq/2.) } ;
  const double eXMY[4] KQED_ALIGN = { exp(-M[0]*xmysq/2.) , exp(-M[1]*xmysq/2.) ,
				      exp(-M[2]*xmysq/2.) , exp(-M[3]*xmysq/2.) } ;

  // look-up table precomputations
  const double Mx[4][4] KQED_ALIGN = { { M[0]*xv[0] , M[1]*xv[0] , M[2]*xv[0] , M[3]*xv[0] } ,
			    { M[0]*xv[1] , M[1]*xv[1] , M[2]*xv[1] , M[3]*xv[1] } ,
			    { M[0]*xv[2] , M[1]*xv[2] , M[2]*xv[2] , M[3]*xv[2] } ,
			    { M[0]*xv[3] , M[1]*xv[3] , M[2]*xv[3] , M[3]*xv[3] } } ;
  const double My[4][4] KQED_ALIGN = { { M[0]*yv[0] , M[1]*yv[0] , M[2]*yv[0] , M[3]*yv[0] } ,
			    { M[0]*yv[1] , M[1]*yv[1] , M[2]*yv[1] , M[3]*yv[1] } ,
			    { M[0]*yv[2] , M[1]*yv[2] , M[2]*yv[2] , M[3]*yv[2] } ,
			    { M[0]*yv[3] , M[1]*yv[3] , M[2]*yv[3] , M[3]*yv[3] } } ;
  const double Mxmy[4][4] KQED_ALIGN = { { M[0]*xmyv[0] , M[1]*xmyv[0] , M[2]*xmyv[0] , M[3]*xmyv[0] } ,
			      { M[0]*xmyv[1] , M[1]*xmyv[1] , M[2]*xmyv[1] , M[3]*xmyv[1] } ,
			      { M[0]*xmyv[2] , M[1]*xmyv[2] , M[2]*xmyv[2] , M[3]*xmyv[2] } ,
			      { M[0]*xmyv[3] , M[1]*xmyv[3] , M[2]*xmyv[3] , M[3]*xmyv[3] } } ;

  // calculate the kernels
  int rhosig , mu , nu , lambda ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {

      // precomputations where I have pulled out the lambda dependence
      // and so these only depend on rhosig and mu indices
      const double Pf2[4] KQED_ALIGN =
	{ yv[0]*Kx0[rhosig][mu][0][0] + yv[1]*Kx0[rhosig][mu][1][0] +
	  yv[2]*Kx0[rhosig][mu][2][0] + yv[3]*Kx0[rhosig][mu][3][0] ,
	  yv[0]*Kx0[rhosig][mu][0][1] + yv[1]*Kx0[rhosig][mu][1][1] +
	  yv[2]*Kx0[rhosig][mu][2][1] + yv[3]*Kx0[rhosig][mu][3][1] ,
	  yv[0]*Kx0[rhosig][mu][0][2] + yv[1]*Kx0[rhosig][mu][1][2] +
	  yv[2]*Kx0[rhosig][mu][2][2] + yv[3]*Kx0[rhosig][mu][3][2] ,
	  yv[0]*Kx0[rhosig][mu][0][3] + yv[1]*Kx0[rhosig][mu][1][3] +
	  yv[2]*Kx0[rhosig][mu][2][3] + yv[3]*Kx0[rhosig][mu][3][3] } ;
      const double Pf3[4] KQED_ALIGN =
	{ yv[0]*K0x[rhosig][0][mu][0] + yv[1]*K0x[rhosig][1][mu][0] +
	  yv[2]*K0x[rhosig][2][mu][0] + yv[3]*K0x[rhosig][3][mu][0] ,
	  yv[0]*K0x[rhosig][0][mu][1] + yv[1]*K0x[rhosig][1][mu][1] +
	  yv[2]*K0x[rhosig][2][mu][1] + yv[3]*K0x[rhosig][3][mu][1] ,
	  yv[0]*K0x[rhosig][0][mu][2] + yv[1]*K0x[rhosig][1][mu][2] +
	  yv[2]*K0x[rhosig][2][mu][2] + yv[3]*K0x[rhosig][3][mu][2] ,
	  yv[0]*K0x[rhosig][0][mu][3] + yv[1]*K0x[rhosig][1][mu][3] +
	  yv[2]*K0x[rhosig][2][mu][3] + yv[3]*K0x[rhosig][3][mu][3] } ;
      const double Pf6[4] KQED_ALIGN =
	{ xmyv[0]*Kx0[rhosig][0][0][mu] + xmyv[1]*Kx0[rhosig][0][1][mu] +
	  xmyv[2]*Kx0[rhosig][0][2][mu] + xmyv[3]*Kx0[rhosig][0][3][mu] ,
	  xmyv[0]*Kx0[rhosig][1][0][mu] + xmyv[1]*Kx0[rhosig][1][1][mu] +
	  xmyv[2]*Kx0[rhosig][1][2][mu] + xmyv[3]*Kx0[rhosig][1][3][mu] ,
	  xmyv[0]*Kx0[rhosig][2][0][mu] + xmyv[1]*Kx0[rhosig][2][1][mu] +
	  xmyv[2]*Kx0[rhosig][2][2][mu] + xmyv[3]*Kx0[rhosig][2][3][mu] ,
	  xmyv[0]*Kx0[rhosig][3][0][mu] + xmyv[1]*Kx0[rhosig][3][1][mu] +
	  xmyv[2]*Kx0[rhosig][3][2][mu] + xmyv[3]*Kx0[rhosig][3][3][mu] } ;
      
      for( nu = 0 ; nu < 4 ; nu++ ) {

	double Pf5 =
	  xv[0]*K0xmy[rhosig][0][nu][mu] +
	  xv[1]*K0xmy[rhosig][1][nu][mu] +
	  xv[2]*K0xmy[rhosig][2][nu][mu] +
	  xv[3]*K0xmy[rhosig][3][nu][mu] ;
	
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {

	  // temporary storage
	  double LMxy[4] KQED_ALIGN , LMyx[4] KQED_ALIGN , LMxxmy[4] KQED_ALIGN ;

	  // inner products go here
	  const double Pf1 =
	    xv[0]*K0y[rhosig][0][nu][lambda] +
	    xv[1]*K0y[rhosig][1][nu][lambda] +
	    xv[2]*K0y[rhosig][2][nu][lambda] +
	    xv[3]*K0y[rhosig][3][nu][lambda] ;
	    
	  const double Pf4 =
	    xv[0]*Ky0[rhosig][nu][0][lambda] +
	    xv[1]*Ky0[rhosig][nu][1][lambda] +
	    xv[2]*Ky0[rhosig][nu][2][lambda] +
	    xv[3]*Ky0[rhosig][nu][3][lambda] ;

	  const double *eXP   = (const double*)eX ;
	  const double *eYP   = (const double*)eY ;
	  const double *eXMYP = (const double*)eXMY ;
	  UNROLL_NM(0);
	  UNROLL_NM(1);
	  UNROLL_NM(2);
	  UNROLL_NM(3); 
	  
	  K[0].L0[rhosig][mu][nu][lambda] = LMxy[0] ;
	  K[1].L0[rhosig][mu][nu][lambda] = LMxxmy[0] ;
	  K[2].L0[rhosig][mu][nu][lambda] = LMxy[0]+LMyx[0]-LMxxmy[0] ;  

	  K[0].L1[rhosig][mu][nu][lambda] = LMxy[1] ;
	  K[1].L1[rhosig][mu][nu][lambda] = LMxxmy[1] ;
	  K[2].L1[rhosig][mu][nu][lambda] = LMxy[1]+LMyx[1]-LMxxmy[1] ;  
	  
	  K[0].L2[rhosig][mu][nu][lambda] = LMxy[2] ;
	  K[1].L2[rhosig][mu][nu][lambda] = LMxxmy[2] ;
	  K[2].L2[rhosig][mu][nu][lambda] = LMxy[2]+LMyx[2]-LMxxmy[2] ;  
	  
	  K[0].L3[rhosig][mu][nu][lambda] = LMxy[3] ;
	  K[1].L3[rhosig][mu][nu][lambda] = LMxxmy[3] ;
	  K[2].L3[rhosig][mu][nu][lambda] = LMxy[3]+LMyx[3]-LMxxmy[3] ;
	}
      }
    }
  }

  return ;
}
