/**
   @file cheby.c
   @brief chebyshev routines
 */
// returns Sum[co[k]*ChebyshevU[k,x],{k,0,nk-1}]
// input: array co[0...(nk-1)], not modified
// assumptions: nk>=1
// implements the Clenshaw Recurrence Formula
// based on Numerical Recipes in C 1992, section 5.5
// U_{n+1}(x) = 2x U_n(x) - U_{n-1}(x)
// i.e.  alpha(n,x)=2*x     beta(n,x)=-1
__device__ KQED_PRIVATE
double
chebUsum( const int nk ,
	  const double x ,
	  const double *co )
{
  double ya=0.0;
  double yb=0.0;
  const double twox = 2.0*x;
  double ytmp ;
  int j;
  for(j=nk-1;j>0;j--) {
    ytmp = yb;
    yb = twox*yb - ya + co[j];
    ya = ytmp;
  }
  return(-ya + twox*yb + co[0]);
}

// Clenshaw for Sum[co[k]*Derivative[0,1][ChebyshevU][k,x],{k,0,nk-1}]
// alf(n,z) = 2z(n+1)/n       beta(n,z) = -(n+2)/n
// y_k = alf(k,x)*y_{k+1} + beta(k+1,x)*y_{k+2} + co[k]
__device__ KQED_PRIVATE
double
dchebUsum( const int nk ,
	   const double x,
	   const double *co)
{
  double ya=0.0;
  double yb=0.0;
  const double twox = 2.0*x;
  double ytmp ;
  int j ;
  for(j=nk-1;j>0;j--) {
    ytmp = yb;
    yb = twox*(j+1)*yb/j - (j+3)*ya/(j+1) + co[j];
    ya = ytmp;
  }
  return(2.0*yb);
}

// Clenshaw for Sum[co[k]*Derivative[0,2][ChebyshevU][k,x],{k,0,nk-1}]
// alf(n,z) = 2z(n+1)/(n-1)       beta(n,z) = -(n+3)/(n-1)
__device__ KQED_PRIVATE
double
ddchebUsum( const int nk,
	    const double x,
	    const double *co)
{
  double ya=0.0;
  double yb=0.0;
  const double twox = 2.0*x;
  double ytmp ;
  int j ;
  for(j=nk-1;j>1;j--) {
    ytmp = yb;
    yb = twox*(j+1)*yb/(j-1) ;
    yb -= (j+4)*ya/j ;
    yb += co[j];
    ya = ytmp;
  }
  return 8*yb ;
}

// Clenshaw for Sum[co[k]*Derivative[0,3][ChebyshevU][k,x],{k,0,nk-1}]
// alf(n,z) = 2z(n+1)/(n-2)       beta(n,z) = -(n+4)/(n-2)
__device__ KQED_PRIVATE
double
dddchebUsum( const int nk ,
	     const double x ,
	     const double *co )
{
  double ya=0.0;
  double yb=0.0;
  const double twox = 2.0*x;
  double ytmp ;
  int j ;
  for(j=nk-1;j>2;j--) {
    ytmp = yb;
    yb = twox*(j+1)*yb/(j-2) - (j+5)*ya/(j-1) + co[j];
    ya = ytmp;
  }
  return 48*yb ;
}
