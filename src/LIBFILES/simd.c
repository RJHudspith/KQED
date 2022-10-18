/**
   @file simd.c
   @brief avx/fma instructions for the hottest parts of the code
   Authors: J. Green and J. Hudspith
 */
#include "KQED.h"

// precompute all this business for Y
void
precompute_INV( struct intprecomp *INVy ,
		const double y ,
		const double y1 ,
		const double y2 ,
		const size_t idx )
{
  INVy -> idx = idx ;
  register const double dy = y1-y2 ;
  register const double ymy1 = y-y1 ;
  register const double ymy2 = y-y2 ;
  register const double ym2sq = (ymy2*ymy2) ;
  register const double ym1sq = (ymy1*ymy1) ;
  INVy -> A = -(ym2sq)*(2*y - 3*y1 + y2) ;
  INVy -> B = (ym1sq)*(2*y + y1 - 3*y2) ;
  INVy -> C1 = (ymy1)*(ym2sq)*dy ;
  INVy -> C2 = (ym1sq)*(ymy2)*dy ;
  INVy -> D = 1./(dy*dy*dy) ;
  INVy -> lA = (y2-y)/(y2-y1) ;
#if (defined HAVE_IMMINTRIN_H) && (defined __AVX__ ) 
  INVy -> y12 = _mm256_setr_pd( y1, y2, y1, y2 ) ;
  INVy -> a  = _mm256_broadcast_sd( &INVy->A  ) ;
  INVy -> b  = _mm256_broadcast_sd( &INVy->B  ) ;
  INVy -> c1 = _mm256_broadcast_sd( &INVy->C1 ) ;
  INVy -> c2 = _mm256_broadcast_sd( &INVy->C2 ) ;
  INVy -> d  = _mm256_broadcast_sd( &INVy->D  ) ;
  INVy -> la = _mm256_broadcast_sd( &INVy->lA  ) ;
#endif
}

#ifdef HAVE_IMMINTRIN_H

#include <immintrin.h>

#ifdef __AVX__

#include "corr_malloc.h" // corr_malloc()
#include "getff-new.h"   // lerp

// look up table space
static __m256d *LUT1 = NULL , *LUT2 = NULL , *LUT3 = NULL ;
static __m256d *LUT4 = NULL , *LUT5 = NULL , *LUT6 = NULL ;

// block inline interpolations
#ifdef __FMA__
#define INBLOCK(A,B)						\
  YMM0 = _mm256_setr_pd( A[0] , A[2] , A[4] , A[6] ) ;		\
  *F = _mm256_mul_pd(YMM0 , INVy.a ) ;				\
  YMM0 = _mm256_setr_pd( A[1] , A[3] , A[5] , A[7] ) ;		\
  *F = _mm256_fmadd_pd( YMM0 , INVy.b , *F ) ;			\
  YMM1 = _mm256_setr_pd( B[0] , B[2] , B[4] , B[6] ) ;		\
  *F = _mm256_fmadd_pd( YMM1 , INVy.c1 , *F ) ;			\
  YMM2 = _mm256_setr_pd( B[1] , B[3] , B[5] , B[7] ) ;		\
  *F = _mm256_fmadd_pd( YMM2 , INVy.c2 , *F ) ;			\
  *F = _mm256_mul_pd( *F , INVy.d ) ;
#else
#define INBLOCK(A,B)						\
  YMM0 = _mm256_setr_pd( A[0] , A[2] , A[4] , A[6] ) ;		\
  *F = _mm256_mul_pd(YMM0 , INVy.a ) ;				\
  YMM0 = _mm256_setr_pd( A[1] , A[3] , A[5] , A[7] ) ;		\
  *F = _mm256_add_pd( *F , _mm256_mul_pd( YMM0 , INVy.b ) ) ;	\
  YMM1 = _mm256_setr_pd( B[0] , B[2] , B[4] , B[6] ) ;		\
  *F = _mm256_add_pd( *F , _mm256_mul_pd( YMM1 , INVy.c1 ) ) ;	\
  YMM2 = _mm256_setr_pd( B[1] , B[3] , B[5] , B[7] ) ;		\
  *F = _mm256_add_pd( *F , _mm256_mul_pd( YMM2 , INVy.c2 ) ) ;	\
  *F = _mm256_mul_pd( *F , INVy.d ) ;
#endif

#ifdef __FMA__

#define CHEB()								\
  ytmp1 = yb1; ytmp2 = yb2;						\
  ya1 = _mm256_sub_pd( ya1 , *co1 ) ;					\
  ya2 = _mm256_sub_pd( ya2 , *co2 ) ;					\
  yb1 = _mm256_fmsub_pd( twox , yb1 , ya1 ) ;				\
  yb2 = _mm256_fmsub_pd( twox , yb2 , ya2 ) ;				\
  co1-- ; co2-- ; ya1 = ytmp1; ya2 = ytmp2;

#define D1CHEB()					\
  ytmp1 = yb1 ; ytmp2 = yb2 ;				\
  ya1 = _mm256_fmsub_pd( *PL2 , ya1 , *co1 ) ;		\
  yb1 = _mm256_mul_pd( yb1 , twox ) ;			\
  ya2 = _mm256_fmsub_pd( *PL2 , ya2 , *co2 ) ;		\
  yb2 = _mm256_mul_pd( yb2 , twox ) ;			\
  yb1 = _mm256_fmsub_pd( *PL1 , yb1 , ya1 ) ;		\
  yb2 = _mm256_fmsub_pd( *PL1 , yb2 , ya2 ) ;		\
  co1 -- ; co2-- ; PL1-- ; PL2-- ;			\
  ya1 = ytmp1 ; ya2 = ytmp2;

#define D2CHEB()					\
  ytmp1 = yb1 ; ytmp2 = yb2 ;				\
  ya1 = _mm256_fmsub_pd( *PL4 , ya1 , *co1 ) ;		\
  yb1 = _mm256_mul_pd( yb1 , twox ) ;			\
  ya2 = _mm256_fmsub_pd( *PL4 , ya2 , *co2 ) ;		\
  yb2 = _mm256_mul_pd( yb2 , twox ) ;			\
  yb1 = _mm256_fmsub_pd( *PL3 , yb1 , ya1 ) ;		\
  yb2 = _mm256_fmsub_pd( *PL3 , yb2 , ya2 ) ;		\
  co1-- ; co2-- ; PL3-- ; PL4-- ;			\
  ya1 = ytmp1 ; ya2 = ytmp2 ;

#define D3CHEB()				\
  ytmp1 = yb1 ; ytmp2 = yb2 ;			\
  ya1 = _mm256_fmsub_pd( *PL6 , ya1 , *co1 ) ;	\
  yb1 = _mm256_mul_pd( yb1 , twox ) ;		\
  ya2 = _mm256_fmsub_pd( *PL6 , ya2 , *co2 ) ;	\
  yb2 = _mm256_mul_pd( yb2 , twox ) ;		\
  yb1 = _mm256_fmsub_pd( *PL5 , yb1 , ya1 ) ;	\
  yb2 = _mm256_fmsub_pd( *PL5 , yb2 , ya2 ) ;	\
  co1-- ; co2-- ; PL5-- ; PL6-- ;		\
  ya1 = ytmp1 ; ya2 = ytmp2 ;			

#else

#define CHEB()						\
  ytmp1 = yb1; ytmp2 = yb2;				\
  yb1 = _mm256_mul_pd( twox , yb1 ) ;			\
  yb2 = _mm256_mul_pd( twox , yb2 ) ;			\
  ya1 = _mm256_sub_pd( ya1 , *co1 ) ;		        \
  ya2 = _mm256_sub_pd( ya2 , *co2 ) ;			\
  yb1 = _mm256_sub_pd( yb1 , ya1 ) ;			\
  yb2 = _mm256_sub_pd( yb2 , ya2 ) ;			\
  co1-- ; co2-- ; ya1 = ytmp1; ya2 = ytmp2;

#define D1CHEB()				\
  ytmp1 = yb1 ; ytmp2 = yb2 ;			\
  yb1 = _mm256_mul_pd( yb1 , twox ) ;		\
  yb2 = _mm256_mul_pd( yb2 , twox ) ;		\
  ya1 = _mm256_mul_pd( ya1 , *PL2 ) ;		\
  ya2 = _mm256_mul_pd( ya2 , *PL2 ) ;		\
  yb1 = _mm256_mul_pd( yb1 , *PL1 ) ;		\
  yb2 = _mm256_mul_pd( yb2 , *PL1 ) ;		\
  ya1 = _mm256_sub_pd( ya1 , *co1 ) ;		\
  ya2 = _mm256_sub_pd( ya2 , *co2 ) ;		\
  yb1 = _mm256_sub_pd( yb1 , ya1 ) ;		\
  yb2 = _mm256_sub_pd( yb2 , ya2 ) ;		\
  co1 -- ; co2-- ; PL1-- ; PL2-- ;		\
  ya1 = ytmp1 ; ya2 = ytmp2;

#define D2CHEB()				\
  ytmp1 = yb1 ; ytmp2 = yb2 ;			\
  yb1 = _mm256_mul_pd( yb1 , twox ) ;		\
  yb2 = _mm256_mul_pd( yb2 , twox ) ;		\
  ya1 = _mm256_mul_pd( ya1 , *PL4 ) ;		\
  ya2 = _mm256_mul_pd( ya2 , *PL4 ) ;		\
  yb1 = _mm256_mul_pd( yb1 , *PL3 ) ;		\
  yb2 = _mm256_mul_pd( yb2 , *PL3 ) ;		\
  ya1 = _mm256_sub_pd( ya1 , *co1 ) ;		\
  ya2 = _mm256_sub_pd( ya2 , *co2 ) ;		\
  yb1 = _mm256_sub_pd( yb1 , ya1 ) ;		\
  yb2 = _mm256_sub_pd( yb2 , ya2 ) ;		\
  co1-- ; co2-- ; PL3-- ; PL4-- ;		\
  ya1 = ytmp1 ; ya2 = ytmp2 ;

#define D3CHEB()				\
  ytmp1 = yb1 ; ytmp2 = yb2 ;			\
  yb1 = _mm256_mul_pd( yb1 , twox ) ;		\
  yb2 = _mm256_mul_pd( yb2 , twox ) ;		\
  ya1 = _mm256_mul_pd( ya1 , *PL6 ) ;		\
  ya2 = _mm256_mul_pd( ya2 , *PL6 ) ;		\
  yb1 = _mm256_mul_pd( yb1 , *PL5 ) ;		\
  yb2 = _mm256_mul_pd( yb2 , *PL5 ) ;		\
  ya1 = _mm256_sub_pd( ya1 , *co1 ) ;		\
  ya2 = _mm256_sub_pd( ya2 , *co2 ) ;		\
  yb1 = _mm256_sub_pd( yb1 , ya1 ) ;		\
  yb2 = _mm256_sub_pd( yb2 , ya2 ) ;		\
  co1-- ; co2-- ; PL5-- ; PL6-- ;		\
  ya1 = ytmp1 ; ya2 = ytmp2 ;

#endif

// turned the call for co1 and co2 into pointer decrements as this
// potentially could allow for loop-unrolling via a look-up-table
static void
chebUsum5( __m256d f[2] ,
	   const int nk ,
	   const double x ,
	   const __m256d *co1,
	   const __m256d *co2 )
{
  const double zero = 0. , TWOx = 2.*x ;
  register __m256d ya1 = _mm256_broadcast_sd( &zero ) ;
  register __m256d yb1 = ya1 , ya2 = ya1 , yb2 = ya1 ;
  register const __m256d twox = _mm256_broadcast_sd( &TWOx ) ;
  register __m256d ytmp1 , ytmp2 ;
  co1 += nk-1 ; co2 += nk-1 ;
#ifdef DUFF
  // loop unroll using duff's device
  register int n = ( (nk) + 7 ) / 8 ;
  switch( (nk)&7 ) {
    case 0 : do { CHEB() ;
    case 7 : CHEB() ;
    case 6 : CHEB() ;
    case 5 : CHEB() ;
    case 4 : CHEB() ;
    case 3 : CHEB() ;
    case 2 : CHEB() ;
    case 1 : CHEB() ;
    } while( --n>0 ) ;
  }
#else
  int i ;
  for( i = 0 ; i < nk ; i++ ) {
    CHEB() ;
  }
#endif
  f[0] = yb1 ;
  f[1] = yb2 ;
}

static void
dchebUsum5( __m256d f[2] ,
	    const int nk ,
	    const double x,
	    const __m256d *co1,
	    const __m256d *co2 )
{
  const double zero = 0. , TWOx = 2.*x , TWO = 2. ;
  const __m256d *PL1 = (const __m256d*)LUT1 ;
  const __m256d *PL2 = (const __m256d*)LUT2 ;
  register __m256d ya1 = _mm256_broadcast_sd( &zero ) ;
  register __m256d yb1 = ya1 , ya2 = ya1 , yb2 = ya2 ;
  register const __m256d twox = _mm256_broadcast_sd( &TWOx ) ;
  co1 += nk-1 ; co2 += nk-1 ;
  PL1 += nk-1 ; PL2 += nk-1 ;
  register __m256d ytmp1 , ytmp2 ;
#ifdef DUFF
  // loop unroll using duff's device
  register int n = ( (nk-1) + 7 ) / 8 ;
  switch( (nk-1)&7 ) {
    case 0 : do { D1CHEB() ;
    case 7 : D1CHEB() ;
    case 6 : D1CHEB() ;
    case 5 : D1CHEB() ;
    case 4 : D1CHEB() ;
    case 3 : D1CHEB() ;
    case 2 : D1CHEB() ;
    case 1 : D1CHEB() ;
    } while( --n>0 ) ;
  }
#else
  int i ;
  for( i = 0 ; i < nk-1 ; i++ ) {
    D1CHEB() ;
  }
#endif
  f[0] = _mm256_mul_pd( yb1 , _mm256_broadcast_sd( &TWO ) ) ;
  f[1] = _mm256_mul_pd( yb2 , _mm256_broadcast_sd( &TWO ) ) ;
}

static void
ddchebUsum5( __m256d f[2] ,
	     const int nk ,
	     const double x,
	     const __m256d *co1,
	     const __m256d *co2 )
{
  const double zero = 0. , TWOx = 2.*x , EIGHT = 8. ;
  const __m256d *PL3 = (const __m256d*)LUT3 ;
  const __m256d *PL4 = (const __m256d*)LUT4 ;
  co1 += nk-1 ; co2 += nk-1 ;
  PL3 += nk-1 ; PL4 += nk-1 ;
  register __m256d ya1 = _mm256_broadcast_sd( &zero ) ;
  register __m256d yb1 = ya1 , ya2 = ya1 , yb2 = ya1 ;
  register const __m256d twox = _mm256_broadcast_sd( &TWOx ) ;
  register __m256d ytmp1 , ytmp2 ;
  // loop unroll using duff's device
#ifdef DUFF
  register int n = ( (nk-2) + 7 ) / 8 ;
  switch( (nk-2)&7 ) {
    case 0 : do { D2CHEB() ;
    case 7 : D2CHEB() ;
    case 6 : D2CHEB() ;
    case 5 : D2CHEB() ;
    case 4 : D2CHEB() ;
    case 3 : D2CHEB() ;
    case 2 : D2CHEB() ;
    case 1 : D2CHEB() ;
    } while( --n>0 ) ;
  }
#else
  int i ;
  for( i = 0 ; i < nk-2 ; i++ ) {
    D2CHEB() ;
  }
#endif
  f[0] = _mm256_mul_pd( yb1 , _mm256_broadcast_sd( &EIGHT ) ) ;
  f[1] = _mm256_mul_pd( yb2 , _mm256_broadcast_sd( &EIGHT ) ) ;
}

static void
dddchebUsum5( __m256d f[2] ,
	      const int nk ,
	      const double x,
	      const __m256d *co1,
	      const __m256d *co2)
{
  const double zero = 0. , TWOx = 2.*x , FEIGHT = 48. ;
  const __m256d *PL5 = (const __m256d*)LUT5 ;
  const __m256d *PL6 = (const __m256d*)LUT6 ;
  co1 += nk-1 ; co2 += nk-1 ;
  PL5 += nk-1 ; PL6 += nk-1 ;
  register __m256d ya1 = _mm256_broadcast_sd( &zero ) ;
  register __m256d yb1 = ya1 , ya2 = ya1 , yb2 = ya1 ;
  register const __m256d twox = _mm256_broadcast_sd( &TWOx ) ;
  register __m256d ytmp1 , ytmp2  ;
  // loop unroll using duff's device
#ifdef DUFF
  register int n = ( (nk-3) + 7 ) / 8 ;
  switch( (nk-3)&7 ) {
    case 0 : do { D3CHEB() ;
    case 7 : D3CHEB() ;
    case 6 : D3CHEB() ;
    case 5 : D3CHEB() ;
    case 4 : D3CHEB() ;
    case 3 : D3CHEB() ;
    case 2 : D3CHEB() ;
    case 1 : D3CHEB() ;
    } while( --n>0 ) ;
  }
#else
  int i ;
  for( i = 0 ; i < nk-3 ; i++ ) {
    D3CHEB() ;
  }
#endif
  f[0] = _mm256_mul_pd( yb1 , _mm256_broadcast_sd( &FEIGHT ) ) ;
  f[1] = _mm256_mul_pd( yb2 , _mm256_broadcast_sd( &FEIGHT ) ) ;
}

// function pointer for cheby stuff
static void (*Func_usm5[4])( __m256d f[2] , const int , const double , const __m256d *co1 , const __m256d *co2 ) = \
{ chebUsum5 , dchebUsum5 , ddchebUsum5 , dddchebUsum5 } ;

static void
getff6( const int *nf, const FFidx nm ,
	const bool calc, const NDCB ndcb, 
	const __m256d y, const double x ,
	const struct AVX_precomps PC[14] ,
	void *FF, void *FFd )
{
  __m256d *ff = (__m256d*)FF ;
  __m256d *ffd = (__m256d*)FFd ;
  
  const int nf_max = (nf[0] > nf[1] ? nf[0] : nf[1]);
  const double ONE = 1.0 ;
  const __m256d one = _mm256_broadcast_sd( &ONE ) ;
  const __m256d y1 = _mm256_div_pd( one , y ) ;

  // enum guarantees these are set but avoid stupid
  // gcc maybe unitialized warning
  __m256d yp ;
  int mshm = 2 ;
  switch(nm) {
  case QG0 : case dxQG0 : yp = one   ; mshm=2 ; break;
  case QG1 : case dxQG1 : yp = one   ; mshm=2 ; break;
  case QG2 : case dxQG2 : case d2xQG2 : yp =  y1   ; mshm=3 ; break;
  case QG3 : case dxQG3 : case d2xQG3 : yp = one   ; mshm=2 ; break;
  case QL4 : case dxQL4 : yp =   y   ; mshm=1 ; break;
  case QL2 : case dxQL2 : yp = _mm256_mul_pd( y1 , y1 ) ; mshm=4 ; break;
  }

  const int offset = nm < dxQG0 ? dxQG0 : QL4  ;
  const int dm = offset + nm > d2xQG3 ? 0 : offset+nm ;

  // these two parameters are simply related to mshm
  int mshp = 2 - mshm ;
  // mm used to be set by the global map Idm, but that was unnecessary
  const int mm = abs( mshp ) ;

  __m256d fval[ nf_max+mm ]  , fvald[ nf_max+mm ];
  __m256d fvaldx[ nf_max+mm ], fvaldxd[ nf_max+mm ];
  
  // call to pow instead of for loop, although mm is expected to be very small
  register __m256d facm = _mm256_mul_pd( y1 , yp ) , facp = yp ;
  facm = _mm256_mul_pd( facm , y1 ) ;

  int j ;
  for( j = 0 ; j < mm ; j++ ) {
    facm = _mm256_mul_pd( facm , y1 );
    facp = _mm256_mul_pd( facp , y ) ;
  }

  // set the remaining parts of fval
  __m256d *Pfval = (__m256d*)fval + mm ;
  __m256d *Pfvald= (__m256d*)fvald+ mm ;
  __m256d *Pfvaldx = (__m256d*)fvaldx + mm ;
  __m256d *Pfvaldxd= (__m256d*)fvaldxd+ mm ;
  
  mshm += mm ;
  mshp += mm ;

  // set MSHM and MSHP here and do everything in double
  __m256d MSHM = _mm256_setr_pd( mshm , mshm , mshm , mshm ) ;
  __m256d MSHP = _mm256_setr_pd( mshp , mshp , mshp , mshp ) ;

  const __m256d *pFm1 = (const __m256d*)PC[nm].Fm ;
  const __m256d *pFm2 = (const __m256d*)PC[dm].Fm ;
  const __m256d *pFp1 = (const __m256d*)PC[nm].Fp ;
  const __m256d *pFp2 = (const __m256d*)PC[dm].Fp ;

#define DO_LOOPS(LOOPB) {						\
    __m256d vm, vp, vmdx, vpdx ;					\
    for( j = 0; j < nf_max; j++ ) {					\
      vm   = _mm256_mul_pd( facm , pFm1[j] ) ;				\
      vp   = _mm256_mul_pd( facp , pFp1[j] ) ;				\
      vmdx = _mm256_mul_pd( facm , pFm2[j] ) ;				\
      vpdx = _mm256_mul_pd( facp , pFp2[j] ) ;				\
      LOOPB ;								\
      Pfval++ ; Pfvald++;						\
      Pfvaldx++ ; Pfvaldxd++;						\
    }									\
  }

  if (calc) {
#ifdef __FMA__
    DO_LOOPS( register __m256d A = _mm256_mul_pd( MSHM , vm ) ;
	      register __m256d B = _mm256_mul_pd( MSHM , vmdx ) ;
	      *Pfval    = _mm256_add_pd( vm , vp ) ;
	      *Pfvaldx  = _mm256_add_pd( vmdx , vpdx ) ;
	      A = _mm256_fmsub_pd( MSHP , vp , A ) ;
	      B = _mm256_fmsub_pd( MSHP , vpdx , B ) ;
	      facm = _mm256_mul_pd( facm , y1 ) ;
	      facp = _mm256_mul_pd( facp , y  ) ;
	      MSHM = _mm256_add_pd( MSHM , one ) ;
	      MSHP = _mm256_add_pd( MSHP , one ) ;
	      *Pfvald   = _mm256_mul_pd( y1 , A ) ;
	      *Pfvaldxd = _mm256_mul_pd( y1 , B ) ; ) ;
#else
    DO_LOOPS( register __m256d A = _mm256_mul_pd( MSHM , vm ) ;
	      register __m256d B = _mm256_mul_pd( MSHP , vp ) ;
	      register __m256d C = _mm256_mul_pd( MSHM , vmdx ) ;
	      register __m256d D = _mm256_mul_pd( MSHP , vpdx ) ;
	      *Pfval = _mm256_add_pd( vm , vp ) ;
	      facm = _mm256_mul_pd( facm , y1 ) ;
	      facp = _mm256_mul_pd( facp , y ) ;
	      A = _mm256_sub_pd( B , A ) ;
	      C = _mm256_sub_pd( D , C ) ;
	      *Pfvaldx  = _mm256_add_pd( vmdx , vpdx ) ;
	      *Pfvald   = _mm256_mul_pd( y1 , A ) ;
	      *Pfvaldxd = _mm256_mul_pd( y1 , C ) ;
	      MSHM = _mm256_add_pd( MSHM , one ) ;
	      MSHP = _mm256_add_pd( MSHP , one ) ; ) ;
#endif
  } else {
    facm = _mm256_mul_pd( facm , y1 ) ;
    facp = _mm256_mul_pd( facp , y1 ) ;
#ifdef __FMA__
    DO_LOOPS( register __m256d A = _mm256_mul_pd( MSHM , vm ) ;
	      register __m256d B = _mm256_mul_pd( MSHM , vmdx ) ;
	      facm = _mm256_mul_pd( facm , y1 ) ;
	      facp = _mm256_mul_pd( facp , y  ) ;
	      MSHM = _mm256_add_pd( MSHM , one ) ;
	      *Pfvald   = _mm256_fmsub_pd( MSHP , vp , A ) ;
	      *Pfvaldxd = _mm256_fmsub_pd( MSHP , vpdx , B  ) ;
	      MSHP = _mm256_add_pd( MSHP , one ) ; ) ;
#else
    DO_LOOPS( register __m256d A = _mm256_mul_pd( MSHM , vm ) ;
	      register __m256d B = _mm256_mul_pd( MSHP , vp ) ;
	      register __m256d C = _mm256_mul_pd( MSHM , vmdx ) ;
	      register __m256d D = _mm256_mul_pd( MSHP , vpdx ) ;
	      facm = _mm256_mul_pd( facm , y1 ) ;
	      facp = _mm256_mul_pd( facp , y  ) ;
	      MSHM = _mm256_add_pd( MSHM , one ) ; 
	      MSHP = _mm256_add_pd( MSHP , one ) ;
	      *Pfvald   = _mm256_sub_pd( B , A ) ;
	      *Pfvaldxd = _mm256_sub_pd( D , C ) ; ) ;
#endif
  }
#undef DO_LOOPS
#undef LOOPA
#undef LOOPC

  if (calc) {
    Func_usm5[mm + ndcb]( ff, nf_max+mm, x, fval, fvaldx ) ;
  }
  Func_usm5[mm + ndcb]( ffd, nf_max+mm, x, fvald, fvaldxd ) ;
}

static void
getff7( const int *nf, const FFidx nm , const NDCB ndcb, 
	const __m256d y, const double x ,
	const struct AVX_precomps PC[14] ,
	void *FF, void *FFd ,
	void *FF2, void *FF2d )
{
  __m256d *ff = (__m256d*)FF ;
  __m256d *ffd = (__m256d*)FFd ;
  __m256d *ff2 = (__m256d*)FF2 ;
  __m256d *ff2d = (__m256d*)FF2d ;
  
  const int nf_max = (nf[0] > nf[1] ? nf[0] : nf[1]);
  const double ONE = 1.0 ;
  const __m256d one = _mm256_broadcast_sd( &ONE ) ;
  const __m256d y1 = _mm256_div_pd( one , y ) ;

  // enum guarantees these are set but avoid stupid
  // gcc maybe unitialized warning
  __m256d yp ;
  int mshm = 2 ;
  switch(nm) {
  case QG0 : case dxQG0 : yp = one   ; mshm=2 ; break;
  case QG1 : case dxQG1 : yp = one   ; mshm=2 ; break;
  case QG2 : case dxQG2 : case d2xQG2 : yp =  y1   ; mshm=3 ; break;
  case QG3 : case dxQG3 : case d2xQG3 : yp = one   ; mshm=2 ; break;
  case QL4 : case dxQL4 : yp =   y   ; mshm=1 ; break;
  case QL2 : case dxQL2 : yp = _mm256_mul_pd( y1 , y1 ) ; mshm=4 ; break;
  }

  const int offset = nm < dxQG0 ? dxQG0 : QL4  ;
  const int dm = offset + nm > d2xQG3 ? 0 : offset+nm ;

  // these two parameters are simply related to mshm
  int mshp = 2 - mshm ;
  // mm used to be set by the global map Idm, but that was unnecessary
  const int mm = abs( mshp ) ;

  __m256d fval[ nf_max+mm ]  , fvald[ nf_max+mm ];
  __m256d fvaldx[ nf_max+mm ], fvaldxd[ nf_max+mm ];
  
  // call to pow instead of for loop, although mm is expected to be very small
  register __m256d facm = _mm256_mul_pd( y1 , yp ) , facp = yp ;
  facm = _mm256_mul_pd( facm , y1 ) ;

  int j ;
  for( j = 0 ; j < mm ; j++ ) {
    facm = _mm256_mul_pd( facm , y1 );
    facp = _mm256_mul_pd( facp , y ) ;
  }

  // set the remaining parts of fval
  __m256d *Pfval    = (__m256d*)fval    + mm ;
  __m256d *Pfvald   = (__m256d*)fvald   + mm ;
  __m256d *Pfvaldx  = (__m256d*)fvaldx  + mm ;
  __m256d *Pfvaldxd = (__m256d*)fvaldxd + mm ;
  
  mshm += mm ;
  mshp += mm ;

  // set MSHM and MSHP here and do everything in double
  __m256d MSHM = _mm256_setr_pd( mshm , mshm , mshm , mshm ) ;
  __m256d MSHP = _mm256_setr_pd( mshp , mshp , mshp , mshp ) ;

  const __m256d *pFm1 = (const __m256d*)PC[nm].Fm ;
  const __m256d *pFm2 = (const __m256d*)PC[dm].Fm ;
  const __m256d *pFp1 = (const __m256d*)PC[nm].Fp ;
  const __m256d *pFp2 = (const __m256d*)PC[dm].Fp ;

#define DO_LOOPS(LOOPB) {						\
    __m256d vm, vp, vmdx, vpdx ;					\
    for( j = 0; j < nf_max; j++ ) {					\
      vm   = _mm256_mul_pd( facm , pFm1[j] ) ;				\
      vp   = _mm256_mul_pd( facp , pFp1[j] ) ;				\
      vmdx = _mm256_mul_pd( facm , pFm2[j] ) ;				\
      vpdx = _mm256_mul_pd( facp , pFp2[j] ) ;				\
      LOOPB ;								\
      Pfval++ ; Pfvald++;						\
      Pfvaldx++ ; Pfvaldxd++;						\
    }									\
  }
  
#ifdef __FMA__
  DO_LOOPS( register __m256d A = _mm256_mul_pd( MSHM , vm ) ;
	    register __m256d B = _mm256_mul_pd( MSHM , vmdx ) ;
	    *Pfval    = _mm256_add_pd( vm , vp ) ;
	    *Pfvaldx  = _mm256_add_pd( vmdx , vpdx ) ;
	    A = _mm256_fmsub_pd( MSHP , vp , A ) ;
	    B = _mm256_fmsub_pd( MSHP , vpdx , B ) ;
	    facm = _mm256_mul_pd( facm , y1 ) ;
	    facp = _mm256_mul_pd( facp , y  ) ;
	    MSHM = _mm256_add_pd( MSHM , one ) ;
	    MSHP = _mm256_add_pd( MSHP , one ) ;
	    *Pfvald   = _mm256_mul_pd( y1 , A ) ;
	    *Pfvaldxd = _mm256_mul_pd( y1 , B ) ; ) ;
#else
  DO_LOOPS( register __m256d A = _mm256_mul_pd( MSHM , vm ) ;
	    register __m256d B = _mm256_mul_pd( MSHP , vp ) ;
	    register __m256d C = _mm256_mul_pd( MSHM , vmdx ) ;
	    register __m256d D = _mm256_mul_pd( MSHP , vpdx ) ;
	    *Pfval = _mm256_add_pd( vm , vp ) ;
	    facm = _mm256_mul_pd( facm , y1 ) ;
	    facp = _mm256_mul_pd( facp , y ) ;
	    A = _mm256_sub_pd( B , A ) ;
	    C = _mm256_sub_pd( D , C ) ;
	    *Pfvaldx  = _mm256_add_pd( vmdx , vpdx ) ;
	    *Pfvald   = _mm256_mul_pd( y1 , A ) ;
	    *Pfvaldxd = _mm256_mul_pd( y1 , C ) ;
	    MSHM = _mm256_add_pd( MSHM , one ) ;
	    MSHP = _mm256_add_pd( MSHP , one ) ; ) ;
#endif
#undef DO_LOOPS
#undef LOOPA
#undef LOOPC

  // always do the derivativw wrt y
  Func_usm5[mm + ndcb]( ff  , nf_max+mm, x , fval  , fvaldx  ) ;
  Func_usm5[mm + ndcb]( ffd , nf_max+mm, x , fvald , fvaldxd ) ;

  // do some more work if we can
  if( (mm + ndcb + 1) < 4 ) {
    Func_usm5[mm + ndcb+1]( ff2  , nf_max+mm, x , fval  , fvaldx  ) ;
    Func_usm5[mm + ndcb+1]( ff2d , nf_max+mm, x , fvald , fvaldxd ) ;
  }

  return ;
}

void
accessv6( const int ix,
	  const FFidx nm,
	  const bool ndy,
	  const NDCB ndcb,
	  const double cb,
	  const struct Grid_coeffs Grid,
	  const struct AVX_precomps *PC ,
	  const struct intprecomp INVy ,
	  double f[4] )
{
  const int nx1 = Grid.nfx[ix];
  const int nx2 = Grid.nfx[ix+1];

  const int nx[2] = { nx1, nx2 };

  double f12[8] __attribute__((aligned(32))) = {0,0,0,0,0,0,0,0} ;
  double g12[8] __attribute__((aligned(32))) = {0,0,0,0,0,0,0,0} ;
  getff6( nx, nm, !ndy, ndcb, INVy.y12, cb, PC, f12, g12 );

  __m256d *F = (__m256d*)f ;
  register __m256d YMM0 , YMM1 , YMM2 ;
  if( ndy ) {
    YMM1 = _mm256_setr_pd( g12[0] , g12[2] , g12[4] , g12[6] ) ;
    YMM2 = _mm256_setr_pd( g12[1] , g12[3] , g12[5] , g12[7] ) ;
    YMM1 = _mm256_sub_pd( YMM1 , YMM2 ) ;
#ifdef __FMA__
    *F = _mm256_fmadd_pd( YMM1 , INVy.la , YMM2 ) ;
#else
    *F = _mm256_add_pd( YMM2 , _mm256_mul_pd( YMM1 , INVy.la ) ) ;
#endif
  } else {
    INBLOCK(f12,g12) ;
  }
  return ;
}

void
accessv7( const int ix,
	  const FFidx nm,
	  const NDCB ndcb,
	  const double cb, 
	  const struct Grid_coeffs Grid,
	  const  struct AVX_precomps *PC ,
	  const struct intprecomp INVy ,
	  double f[12] )
{
  const int nx1 = Grid.nfx[ix];
  const int nx2 = Grid.nfx[ix+1];

  const int nx[2] = { nx1, nx2 };

  double f12[8] __attribute__((aligned(32))) = {0,0,0,0,0,0,0,0} ;
  double g12[8] __attribute__((aligned(32))) = {0,0,0,0,0,0,0,0} ;
  double f212[8] __attribute__((aligned(32))) = {0,0,0,0,0,0,0,0} ;
  double g212[8] __attribute__((aligned(32))) = {0,0,0,0,0,0,0,0} ;

  getff7( nx, nm, ndcb, INVy.y12, cb, PC, f12, g12 , f212 , g212 );

  __m256d *F = (__m256d*)f ;
  register __m256d YMM0 , YMM1 , YMM2 ;
  
  INBLOCK(f12,g12) ; F++ ;
  YMM1 = _mm256_sub_pd( YMM1 , YMM2 ) ;
#ifdef __FMA__
  *F = _mm256_fmadd_pd( YMM1 , INVy.la , YMM2 ) ;
#else
  *F = _mm256_add_pd( YMM2 , _mm256_mul_pd( YMM1 , INVy.la ) ) ;
#endif
  F++ ;
  INBLOCK(f212,g212) ;
  return ;
}

// free our allocated look-up-tables (LUTs)
void
free_LUTs( void )
{
  if( LUT1 != NULL ) {
    free( LUT1 ) ;
  }
  if( LUT2 != NULL ) {
    free( LUT2 ) ;
  }
  if( LUT3 != NULL ) {
    free( LUT3 ) ;
  }
  if( LUT4 != NULL ) {
    free( LUT4 ) ;
  }
  if( LUT5 != NULL ) {
    free( LUT5 ) ;
  }
  if( LUT6 != NULL ) {
    free( LUT6 ) ;
  }
  return ;
}

// compute the first two LUTs
int
init_LUT12( const size_t max_LUT )
{
  if( corr_malloc( (void**)&LUT1 , 32 , max_LUT*sizeof(__m256d) ) != 0 ||
      corr_malloc( (void**)&LUT2 , 32 , max_LUT*sizeof(__m256d) ) != 0 ) {
    fprintf( stderr , "[INIT] LUT12 alloc failed\n" ) ;
    return 1 ;
  }
  double tmp = 0 ;
  LUT1[0] = _mm256_broadcast_sd( &tmp ) ; 
  LUT2[0] = _mm256_broadcast_sd( &tmp ) ;
  size_t i ;
  for( i = 1 ; i < max_LUT ; i++ ) {
    tmp = (i+1)/(double)(i) ;
    LUT1[i] = _mm256_broadcast_sd( &tmp ) ; 
    tmp = (i+3)/(double)(i+1) ;
    LUT2[i] = _mm256_broadcast_sd( &tmp ) ; 
  }
  return 0 ;
}

// compute the next two LUTs
int
init_LUT34( const size_t max_LUT )
{
  if( corr_malloc( (void**)&LUT3 , 32 , max_LUT*sizeof(__m256d) ) != 0 ||
      corr_malloc( (void**)&LUT4 , 32 , max_LUT*sizeof(__m256d) ) != 0 ) {
    fprintf( stderr , "[INIT] LUT34 alloc failed\n" ) ;
    return 1 ;
  }
  double tmp = 0 ;
  LUT3[0] = _mm256_broadcast_sd( &tmp ) ; 
  LUT4[0] = _mm256_broadcast_sd( &tmp ) ; 
  LUT3[1] = _mm256_broadcast_sd( &tmp ) ; 
  LUT4[1] = _mm256_broadcast_sd( &tmp ) ; 
  size_t i ;
  for( i = 2 ; i < max_LUT ; i++ ) {
    tmp = (i+1)/(double)(i-1) ;
    LUT3[i] = _mm256_broadcast_sd( &tmp ) ;
    tmp = (i+4)/(double)(i) ;
    LUT4[i] = _mm256_broadcast_sd( &tmp ) ; 
  }
  return 0 ;
}

// and the final two
int
init_LUT56( const size_t max_LUT )
{
  if( corr_malloc( (void**)&LUT5 , 32 , max_LUT*sizeof(__m256d) ) != 0 ||
      corr_malloc( (void**)&LUT6 , 32 , max_LUT*sizeof(__m256d) ) != 0 ) {
    fprintf( stderr , "[INIT] LUT56 alloc failed\n" ) ;
    return 1 ;
  }
  double tmp = 0 ;
  LUT5[0] = _mm256_broadcast_sd( &tmp ) ; 
  LUT6[0] = _mm256_broadcast_sd( &tmp ) ; 
  LUT5[1] = _mm256_broadcast_sd( &tmp ) ; 
  LUT6[1] = _mm256_broadcast_sd( &tmp ) ; 
  LUT5[2] = _mm256_broadcast_sd( &tmp ) ; 
  LUT6[2] = _mm256_broadcast_sd( &tmp ) ; 
  size_t i ;
  for( i = 3 ; i < max_LUT ; i++ ) {
    tmp = (i+1)/(double)(i-2) ;
    LUT5[i] = _mm256_broadcast_sd( &tmp ) ;
    tmp = (i+5)/(double)(i-1) ;
    LUT6[i] = _mm256_broadcast_sd( &tmp ) ; 
  }
  return 0 ;
}

#endif
#endif // HAVE_IMMINTRIN_H
