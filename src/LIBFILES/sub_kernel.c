/**
   @file sub_kernel.c
   @brief computes the kernel needed for the subleading contributions (3+1 and 2+1+1)

   ( dropping rho,sigma arguments)
   
   A = L_{\mu\nu\lambda}(x,y)
   B = L_{\nu\lambda\mu}(x,x-y)
   C = L_{\mu\nu\lambda}(x,y) + L_{\nu\mu\lambda}(y,x) - L_{\nu\lambda\mu}(x,x-y)
*/
#include "KQED.h"

#define UNROLL_NM_SUB(Nm)							\
  LMxy[Nm] = Kxy[rhosig][mu][nu][lambda]				\
    -*eXP*( K0y[rhosig][mu][nu][lambda] - Mx[mu][Nm]*Pf1 )		\
    -*eYP*( Kx0[rhosig][mu][nu][lambda] - My[nu][Nm]*Pf2[lambda] )  ;	\
  LMyx[Nm] = Kyx[rhosig][nu][mu][lambda]				\
    -*eYP*( K0x[rhosig][nu][mu][lambda]-My[nu][Nm]*Pf3[lambda] )	\
    -*eXP*( Ky0[rhosig][nu][mu][lambda]-Mx[mu][Nm]*Pf4 )  ;		\
  LMxxmy[Nm] = Kxxmy[rhosig][nu][lambda][mu]				\
    -*eXP*( K0xmy[rhosig][nu][lambda][mu]-Mx[nu][Nm]*Pf5[lambda] )		\
    -*eXMYP*( Kx0[rhosig][nu][lambda][mu]-Mxmy[lambda][Nm]*Pf6 )  ; \
  eXP++ ; eYP++ ; eXMYP++ ;

// compute 4 L2 kernels for the subleading diagrams with the M-factor included
void
compute_sub_kernelsM_L2( const double M[4] ,
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
      const double Pf5[4] KQED_ALIGN =
	{ xv[0]*K0xmy[rhosig][0][0][mu] + xv[1]*K0xmy[rhosig][1][0][mu] +
	  xv[2]*K0xmy[rhosig][2][0][mu] + xv[3]*K0xmy[rhosig][3][0][mu] ,
	  xv[0]*K0xmy[rhosig][0][1][mu] + xv[1]*K0xmy[rhosig][1][1][mu] +
	  xv[2]*K0xmy[rhosig][2][1][mu] + xv[3]*K0xmy[rhosig][3][1][mu] ,
	  xv[0]*K0xmy[rhosig][0][2][mu] + xv[1]*K0xmy[rhosig][1][2][mu] +
	  xv[2]*K0xmy[rhosig][2][2][mu] + xv[3]*K0xmy[rhosig][3][2][mu] ,
	  xv[0]*K0xmy[rhosig][0][3][mu] + xv[1]*K0xmy[rhosig][1][3][mu] +
	  xv[2]*K0xmy[rhosig][2][3][mu] + xv[3]*K0xmy[rhosig][3][3][mu] } ;
      
      for( nu = 0 ; nu < 4 ; nu++ ) {
	const double Pf6 =
	  xmyv[0]*Kx0[rhosig][nu][0][mu] +
	  xmyv[1]*Kx0[rhosig][nu][1][mu] +
	  xmyv[2]*Kx0[rhosig][nu][2][mu] +
	  xmyv[3]*Kx0[rhosig][nu][3][mu] ;
	
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
	  UNROLL_NM_SUB(0);
	  UNROLL_NM_SUB(1);
	  UNROLL_NM_SUB(2);
	  UNROLL_NM_SUB(3); 
	  
	  K[0].L0[rhosig][mu][nu][lambda] = LMxy[0] ;
	  K[1].L0[rhosig][mu][nu][lambda] = LMxxmy[0] ;
	  K[2].L0[rhosig][mu][nu][lambda] = LMxy[0] + LMyx[0] - LMxxmy[0] ;  

	  K[0].L1[rhosig][mu][nu][lambda] = LMxy[1] ;
	  K[1].L1[rhosig][mu][nu][lambda] = LMxxmy[1] ;
	  K[2].L1[rhosig][mu][nu][lambda] = LMxy[1] + LMyx[1] - LMxxmy[1] ;  
	  
	  K[0].L2[rhosig][mu][nu][lambda] = LMxy[2] ;
	  K[1].L2[rhosig][mu][nu][lambda] = LMxxmy[2] ;
	  K[2].L2[rhosig][mu][nu][lambda] = LMxy[2] + LMyx[2] - LMxxmy[2] ;  
	  
	  K[0].L3[rhosig][mu][nu][lambda] = LMxy[3] ;
	  K[1].L3[rhosig][mu][nu][lambda] = LMxxmy[3] ;
	  K[2].L3[rhosig][mu][nu][lambda] = LMxy[3] + LMyx[3] - LMxxmy[3] ;
	}
      }
    }
  }

  return ;
}
