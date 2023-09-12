#include "KQED.h"      // definitions and enums

#include "cheby.h"     // chebUsum and alike
#include "getff-new.h" // alphabetising

// interpolation function
// dy should be set to y1-y2
__device__
static inline double
interpol3( const double y ,
	   const double y1 ,
	   const double y2 ,
	   const double dy ,
	   const double f1 ,
	   const double f2 ,
	   const double g1 ,
	   const double g2 )
{
  return ((-f1* (y - y2)*(y - y2)*(2* y - 3* y1 + y2) +			\
	   (y - y1)* (f2* (y - y1)*(2* y + y1 - 3* y2)			\
		      + (y - y2)*dy* (g1* y + g2* y - g2* y1 - g1* y2)))/(dy*dy*dy)) ;
}

// precompute all this business for x or y
__device__ KQED_PRIVATE
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
}

// precompute all this business for x
__device__ KQED_PRIVATE
void
precompute_INVx( struct intprecomp *INVx ,
    const double y ,
    const double y1 ,
    const double y2 ,
    const size_t idx )
{
  precompute_INV(INVx, y, y1, y2, idx);
}

__device__
static inline double
interpol4( const struct intprecomp INVx ,
	   const double f1 ,
	   const double f2 ,
	   const double g1 ,
	   const double g2 )
{
  return (f1*INVx.A + f2*INVx.B + g1*INVx.C1 + g2*INVx.C2)*INVx.D ; 
}

// function pointer for cheby stuff
// static double (*Func_usm[4])( const int , const double , const double *) =
// { chebUsum , dchebUsum , ddchebUsum , dddchebUsum } ;
__device__
static double Func_usm( const int i,
    const int n, const double x, const double *f ) {
  switch(i) {
    case 0: return chebUsum(n, x, f);
    case 1: return dchebUsum(n, x, f);
    case 2: return ddchebUsum(n, x, f);
    case 3: return dddchebUsum(n, x, f);
    default: return 0.0;
  }
}

// returns the form factor, given the coefficients fm[0..(nf-1)] and fp[0..(nf-1)]
// e.g. fm = alpha^{(3)}_{m-}  and fp = alpha^{(3)}_{m+}
// nf = length of vectors fm[nm][ix/iy] and similarly for fp
// nm = index of the sum sigma that appeared in the integrand == outer index of ff
// ndy and ndcb = # derivatives with respect to y and cos(beta) respectively
// x is not used in this function
__device__
static void
getff2( double res[2] ,
	const int nf, const FFidx nm ,
	const bool ndy , const NDCB ndcb, 
	const double y, const double x ,
	const float *fm, const float *fp )
{
  const double y1 = 1.0/y;

  // enum guarantees these are set but avoid stupid
  // gcc maybe unitialized warning
  double yp = 1 ;
  int mshm = 2 ;
  switch(nm) {
  case QG0    : case dxQG0 : yp = 1.0   ; mshm=2 ; break;
  case QG1    : case dxQG1 : yp = 1.0   ; mshm=2 ; break;
  case QG2    : case dxQG2 : case d2xQG2: yp =  y1 ; mshm=3 ; break;
  case QG3    : case dxQG3 : case d2xQG3: yp = 1.0 ; mshm=2 ; break;
  case QL4    : case dxQL4 : yp =   y   ; mshm=1 ; break;
  case QL2    : case dxQL2 : yp = y1*y1 ; mshm=4 ; break;
  }

  // these two parameters are simply related to mshm
  int mshp = 2 - mshm ;
  // mm used to be set by the global map Idm, but that was unnecessary
  const int mm = abs( mshp ) ;

  // set fval to zero, actually not really needed
  double fval[128], fvalD[128];
  // TODO: do something more reasonable if too long
  if (nf+mm >= 128) return;
  memset( fval , 0 , (nf+mm)*sizeof( double ) ) ;
  memset( fvalD , 0 , (nf+mm)*sizeof( double ) ) ;
  double facm = y1*y1*yp ;
  double facp = yp ;

  int j ;
  for( j = 0 ; j < mm ; j++ ) {
    facm *= y1 ;
    facp *= y ;
  }

  double *Pfval = (double*)fval + mm ;
  double *PfvalD = (double*)fvalD + mm ;

  mshm += mm ; mshp += mm ;
  for(j=0;j<nf;j++) {
    *PfvalD = y1*mshp*facp*(*fp) - y1*mshm*facm*(*fm) ;
    *Pfval = facm*(*fm) + facp*(*fp) ;
    facm *= y1; facp *= y;    
    Pfval++ ; PfvalD++ ; fm++ ; fp++ ; mshm++ ; mshp++ ;
  }

  res[1] = Func_usm(mm + ndcb, nf+mm, x , fvalD ) ;
  res[0] = ndy? res[1] : Func_usm(mm + ndcb, nf+mm, x , fval ) ;
}

// case where you have read in the weight functions upon initialization
// interpolates the form factor[nm] to the target point y using the grid
__device__ KQED_PRIVATE
double
accessv( const bool flag_hy, const bool use_y_derivs,
	 const int ix, const int iy,
	 const FFidx nm, const bool ndy, const NDCB ndcb,
	 const double cb, const double y, 
	 const struct Grid_coeffs Grid )
{
  const int nx = Grid.nfx[ix] ;
  const double y1 = Grid.YY[iy];
  double res1[2] = {0.,0.} ;
    
  // occasionally ndy gets set to 1, if use_y_derivs is set we always do the deriv
  getff2( res1 , nx, nm, ndy, ndcb, y1, cb,
	  getFfm(&Grid, nm, ix, iy) , getFfp(&Grid, nm, ix, iy) ) ;
  
  // if we are not at the upper limit of Y we can use info from the next point
  if(!flag_hy) {
    const int iy2 = iy+1;
    const double y2 = Grid.YY[iy2];
    double res2[2] = {0.,0.} ;

    getff2( res2 , nx, nm, ndy, ndcb, y2, cb,
	    getFfm(&Grid, nm, ix, iy2) ,
	    getFfp(&Grid, nm, ix, iy2) ) ;
    
    if( use_y_derivs ) {
      return interpol3( y, y1, y2, y1-y2,
			res1[0], res2[0], res1[1], res2[1] ) ;      
    } else {
      // value interpolated to target y, at x=x1[0];
      return lerp( (y2-y)/(y2-y1) , res1[0] , res2[0] ) ;
    }
  }
  return res1[0] ;
}

// returns the lower index that bounds "target"
// e.g arr[lo] < target < arr[lo+1]
// (assumes a monotonically increasing arr)
__device__
static int
bsrch( const double *arr , const double target ,
  const int lo , const int hi )
{
  // when hi == lo we are done
  if( ( hi - lo ) < 2 ) return lo ;
  const int mid = ( hi + lo )/2 ;
  if( arr[mid] > target ) {
    return bsrch( arr , target , lo , mid ) ;
  } else {
    return bsrch( arr , target , mid , hi ) ;
  }
}

// linear search variant targeted at avoiding warp divergence and keeping
// regular memory access on CUDA.
__device__
static int
lsrch( const double *arr, const double target,
    const int lo, const int hi ) {
  int index = lo;
  for (int i = lo; i < hi; ++i) {
    index = (arr[i] <= target && target <= arr[i+1]) ? i : index;
  }
  return index;
}

// returns the lower index that bounds "target"
// e.g arr[lo] < target < arr[lo+1]
// (assumes a monotonically increasing arr)
__device__ KQED_PRIVATE
int
find_ind(const double *arr, const double target,
    const int lo, const int hi) {
  // NOTE: this is not well-suited for GPU execution
  // return bsrch(arr, target, lo, hi);
  // Better option:
  return lsrch(arr, target, lo, hi);
}

// extract the form factor
__device__ KQED_PRIVATE
double
extractff( const FFidx nm, const bool ndy, const NDCB ndcb,
	   const struct invariants Inv , const struct Grid_coeffs Grid )
{  
  const bool flag_hx = ( Inv.x >= Grid.XX[ Grid.nstpx-1 ] ) ;
  const bool flag_hy = ( Inv.y >= Grid.YY[ Grid.nstpy-1 ] ) ;
 
  const bool use_x_derivs = (nm<dxQG0 || nm==dxQG2 || nm==dxQG3) ;

  const int ix = Inv.INVx.idx , iy = Inv.INVy.idx ; 
  const bool use_y_derivs = (ndy==false) ;
  const double f1iy = accessv( flag_hy, use_y_derivs, ix, iy, nm,
			       ndy, ndcb, Inv.cb, Inv.y, Grid ) ;

  if(!flag_hx) {
    const double f2iy = accessv( flag_hy, use_y_derivs, ix+1, iy, nm,
				 ndy, ndcb, Inv.cb, Inv.y , Grid ) ;
    if(use_x_derivs) {
      const int offset = nm < dxQG0 ? dxQG0 : QL4 ;
      // These enum additions are pretty sketchy...
      const double g1iy = accessv( flag_hy, use_y_derivs, ix, iy, (FFidx)(nm+offset),
				   ndy, ndcb, Inv.cb, Inv.y, Grid ) ;
      const double g2iy = accessv( flag_hy, use_y_derivs, ix+1, iy, (FFidx)(nm+offset),
				   ndy, ndcb, Inv.cb, Inv.y, Grid ) ;

      return interpol4( Inv.INVx, f1iy, f2iy, g1iy, g2iy ) ;
    } else {
      // lerpity lerp
      return lerp( Inv.INVx.lA , f1iy , f2iy ) ;
    }
  }
  return f1iy ;
}

// extract the form factor
__device__ KQED_PRIVATE
void
extractff2( const FFidx nm,
	    const NDCB ndcb,
	    const struct invariants Inv ,
	    const struct Grid_coeffs Grid ,
	    double F[4] )
{  
  // derivative map
  const NDCB ndcb2 = (NDCB)(ndcb+1 < 5 ? ndcb+1 : ndcb) ;
  // map for the x-derivative, incomplete as some derivatives aren't used
  const FFidx dxmap[14] = { dxQG0  , dxQG1 , dxQG2  , dxQG3  , dxQL4 , dxQL2 ,
			  dxQG0  , dxQG1 , d2xQG2 , d2xQG3 , dxQL4 , dxQL2 ,
			  d2xQG2 , d2xQG3 } ;
  // derivative wrt dcb
  F[0] = extractff( nm, false , ndcb2 , Inv , Grid ) ;
  // derivative wrt x
  F[1] = extractff( dxmap[nm], false , ndcb  , Inv , Grid ) ;
  // derivative wrt y
  F[2] = extractff( nm, true  , ndcb  , Inv , Grid ) ;
  // no derivative
  F[3] = extractff( nm, false , ndcb  , Inv , Grid ) ;
  return ;
}

// Initialises the invariants used in chnr_*
__device__ KQED_PRIVATE
struct invariants
set_invariants( const double xv[4] ,
		const double yv[4] ,
		const struct Grid_coeffs Grid )
{
  struct invariants Inv ;
  const double EPSIN = 1E-7 ;
  
  Inv.xsq = SCALPROD( xv , xv ) ;
  Inv.xsq = Inv.xsq < EPSIN ? EPSIN : Inv.xsq ;
  Inv.x = sqrt( fabs( Inv.xsq ) ) ;

  // y needs a little fudge factor as it is badly behaved for very small y
  Inv.ysq = SCALPROD( yv , yv ) ;
  Inv.ysq = Inv.ysq < EPSIN ? EPSIN : Inv.ysq ;
  Inv.y = sqrt( fabs( Inv.ysq ) ) ;

  Inv.xdoty = SCALPROD( xv , yv ) ;
  Inv.cb = Inv.xdoty/( Inv.x * Inv.y ) ;

  Inv.cborig = Inv.cb ;
  Inv.yorig = Inv.y ;
  Inv.xmysq = ( Inv.xsq + Inv.ysq )*(1.00000000000001)-2.0*Inv.xdoty ;
  Inv.xmy = sqrt( fabs( Inv.xmysq ) ) ;

  Inv.flag2 = false ;
  const double rx = ( Inv.x > Inv.xmy ? Inv.xmy/Inv.x : Inv.x/Inv.xmy);
  const double rxy = ( Inv.x > Inv.y ? Inv.y/Inv.x : Inv.x/Inv.y);

  if( rx < rxy
      && Inv.xmy < Grid.YY[ Grid.nstpy-1] 
      && fabs(Inv.xmysq) > 1E-28 ) {
    Inv.cb = (Inv.x-Inv.y*Inv.cb)/Inv.xmy ;
    Inv.y = Inv.xmy ;    
    Inv.flag2 = true ;
  }

  // setup INVx
  const size_t ix1 = (size_t)find_ind( Grid.XX , Inv.x , 0 , Grid.nstpx ) ;
  size_t ix2 = ix1+1 ;
  if( ix2 >= (size_t)Grid.nstpx ) {
    ix2 = ix1 ;
  }
  precompute_INVx( &Inv.INVx, Inv.x, Grid.XX[ix1], Grid.XX[ix2] , ix1 ) ;
  
  // setup InvY
  const size_t iy1 = (size_t)find_ind( Grid.YY , Inv.y , 0 , Grid.nstpy ) ;  
  // y edge case
  size_t iy2 = iy1 + 1 ;
  if( iy2 >= (size_t)Grid.nstpy ) {
    iy2 = iy1 ;
  }
  precompute_INV( &Inv.INVy , Inv.y , Grid.YY[iy1] , Grid.YY[iy2] , iy1 ) ;
  
  return Inv ;
}
