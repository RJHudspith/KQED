#ifndef KERNELS_H
#define KERNELS_H

// from linear idx in 0->384 return a value where mu and lambda indices are swapped
__device__
static inline size_t
i_to_mulam( const size_t idx )
{
  const size_t l[4] = { idx/64 , (idx/16)&3 , (idx/4)&3 , idx&3 } ;
  return l[1] + 4*(l[2]+4*(l[3]+4*l[0]) ) ;
}


__device__
void
QED_kernel_L0( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

__device__
void
QED_kernel_L1( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

__device__
void
QED_kernel_L2( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

__device__
void
QED_Mkernel_L2( const double M ,
		const double xv[4] ,
		const double yv[4] ,
		const struct QED_kernel_temps t ,
		double kerv[6][4][4][4] ) ;

__device__
void
QED_kernel_L3( const double xv[4] ,
	       const double yv[4] ,
	       const struct QED_kernel_temps t ,
	       double kerv[6][4][4][4] ) ;

#endif
