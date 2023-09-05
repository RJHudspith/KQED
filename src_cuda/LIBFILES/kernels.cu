/**
   @file kernels.c
   @brief kernel codes
 */
#include "KQED.h"

#include "kernels.h"        // alphabetising
#include "QED_kernel.h"     // kernelQED()
#include "QED_kernel_xy0.h" // kernelQED_xoryeq0()
#include "Tabd.h"

// check that a 4-vector is zero
__device__
static inline bool
is_zero( const double x[4] )
{
  return (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]) < 1E-28 ;
}

// check if x vector is equal to y vector
__device__
static inline bool
x_is_y( const double x[4] ,
	const double y[4] )
{
  const double xmy[4] = {x[0]-y[0],x[1]-y[1],x[2]-y[2],x[3]-y[3]} ;
  return is_zero( xmy ) ;
}

// from linear idx in 0->384 return a value where mu and lambda indices are swapped
__device__
static inline size_t
i_to_mulam( const size_t idx )
{
  const size_t l[4] = { idx/64 , (idx/16)&3 , (idx/4)&3 , idx&3 } ;
  return l[1] + 4*(l[2]+4*(l[3]+4*l[0]) ) ;
}

// sets L -> -L_{mu<->lambda}
__device__
static void
kernel_mulam_minus( double *kerv )
{
  double ktmp[384] ;
  memcpy( ktmp , kerv , 384*sizeof(double) ) ;
  size_t i ;
  for( i = 0 ; i < 384 ; i++ ) {
    *kerv = -ktmp[ i_to_mulam(i) ] ;
    kerv++ ;
  }
  return ;
}

// computes kerv += kt * S index-by-index
__device__
static void
atomic_kernel_Saxpy( double *kerv ,
		     const double S ,
		     const double *kt )
{
  int i ;
  // although not strictly necessary is a nice proof that we
  // are correctly aligned to the boundary
  for( i = 0 ; i < 384 ; i++ ) {
    *kerv += *kt * S ;
    kerv++ ; kt++ ;
  }
  return ;
}

// L0 is the standard version of the kernel
__device__
void
QED_kernel_L0( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] )
{
  const bool x_is_zero  = is_zero( xv ) ;
  const bool y_is_zero  = is_zero( yv ) ;

  // set kernel to zero
  memset( kerv , 0 , 384*sizeof( double ) ) ;
  
  // do the logic
  if( x_is_zero && y_is_zero ) {
    return ;
  }
  if( x_is_zero ) {
    kernelQED_xoryeq0( yv , t , kerv , Tabd_xeq0 ) ;
    return ;
  }
  if( y_is_zero ) {
    kernelQED_xoryeq0( xv , t , kerv , Tabd_yeq0 ) ;
    return ;
  }
  if( x_is_y( xv , yv ) ) {
    kernelQED_xoryeq0( xv , t , kerv , Tabd_yeq0 ) ;
    kernel_mulam_minus( (double*)kerv ) ;
    return ;
  }
  
  // if neither are zero or small we do the standard one
  kernelQED( xv , yv , t , kerv ) ;
  
  return ;
}

// L1 = L0(x,y) - L0(x,x)/2 - L0(y,y)/2
__device__
void
QED_kernel_L1( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] )
{
  double kt[6][4][4][4] KQED_ALIGN ;
  QED_kernel_L0( xv , yv , t , kerv ) ;

  QED_kernel_L0( xv , xv , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , -0.5 , (const double*)kt ) ; 

  QED_kernel_L0( yv , yv , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , -0.5 , (const double*)kt ) ;
  
  return ;
}

// L2 = L0(x,y) - L0(0,y) - L0(x,0)
__device__
void
QED_kernel_L2( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] )
{
  double kt[6][4][4][4] KQED_ALIGN ;
  const double zero[4] = { 0 , 0 , 0 , 0 } ;
  QED_kernel_L0( xv , yv , t , kerv ) ;
  
  QED_kernel_L0( zero , yv , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , -1 , (const double*)kt ) ;
  
  QED_kernel_L0( xv , zero , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , -1 , (const double*)kt ) ;
  
  return ;
}

// L3 = L0(x,y) - L0(x,x) + L0(0,x) - L0(0,y)
__device__
void
QED_kernel_L3( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] )
{
  double kt[6][4][4][4] KQED_ALIGN ;
  const double zero[4] = { 0 , 0 , 0 , 0 } ;
  QED_kernel_L0( xv , yv , t , kerv ) ;

  QED_kernel_L0( xv , xv , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , -1 , (const double*)kt ) ;

  QED_kernel_L0( zero , xv , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , +1 , (const double*)kt ) ;

  QED_kernel_L0( zero , yv , t , kt ) ;
  atomic_kernel_Saxpy( (double*)kerv , -1 , (const double*)kt ) ;
  
  return ;
}

// L2 = L0(x,y) - d_\mu(x_\alpha exp^{-Mx^2/2}) [L0(0,y) - L0(x,0)]
__device__
void
QED_Mkernel_L2( const double M ,
		const double xv[4] ,
		const double yv[4] ,
		const struct QED_kernel_temps t ,
		double kerv[6][4][4][4] )
{   
  double kt[6][4][4][4] KQED_ALIGN ;
  const double zero[4] = { 0 , 0 , 0 , 0 } ;
  QED_kernel_L0( xv , yv , t , kerv ) ;

  int rhosig , mu , nu , lambda ;
  const double gaussX = exp( -M*(xv[0]*xv[0]+xv[1]*xv[1]+
				 xv[2]*xv[2]+xv[3]*xv[3] )/2. ) ;

  QED_kernel_L0( zero , yv , t , kt ) ;
  
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  kerv[rhosig][mu][nu][lambda] -=
	    gaussX*( kt[rhosig][mu][nu][lambda] -
		     M*xv[mu]*( xv[0]*kt[rhosig][0][nu][lambda] +
				xv[1]*kt[rhosig][1][nu][lambda] +
				xv[2]*kt[rhosig][2][nu][lambda] +
				xv[3]*kt[rhosig][3][nu][lambda] ) ) ;
	}
      }
    }
  }

  const double gaussY = exp( -M*(yv[0]*yv[0]+yv[1]*yv[1]+
				 yv[2]*yv[2]+yv[3]*yv[3] )/2. ) ;
  QED_kernel_L0( xv , zero , t , kt ) ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	  kerv[rhosig][mu][nu][lambda] -=
	    gaussY*( kt[rhosig][mu][nu][lambda] -
		     M*yv[nu]*( yv[0]*kt[rhosig][mu][0][lambda] +
				yv[1]*kt[rhosig][mu][1][lambda] +
				yv[2]*kt[rhosig][mu][2][lambda] +
				yv[3]*kt[rhosig][mu][3][lambda] ) ) ;
	}
      }
    }
  }
  
  return ;
}
