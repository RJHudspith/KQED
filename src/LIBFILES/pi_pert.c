/**
   @file pi_pert.c
   @brief various perturbative results for pi
 */
#include "KQED.h"    // definitions and what have you 

#include "pi_pert.h" // alphabetising

#define pi M_PI
#define ZETA3 (1.2020569031595942854)

// norm of a vector
static double
v4norm( const double *x )
{
  return sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] ) ;
}

// kronecker delta
static inline int
kro( const int mu , const int nu )
{
  return (mu==nu) ;
}

// lexicographical order for G8
static inline int
lexi_G8( const int a , const int b , const int c , const int d ,
	 const int e , const int f , const int g , const int h )
{
  return h+4*(g+4*(f+4*(e+4*(d+4*(c+4*(b+4*a)))))) ;
}

// lexicographical order for G8
static inline int
lexi_G6( const int c , const int d ,
	 const int e , const int f , const int g , const int h )
{
  return h+4*(g+4*(f+4*(e+4*(d+4*c)))) ;
}

// lexicographical order for G4
static inline int
lexi_G4( const int e , const int f , const int g , const int h )
{
  return h+4*(g+4*(f+4*e)) ;
}

// returns contribution of the lbl contribution due to a heavy lepton; mh = (lepton mass)/(muon mass)
// including the (alpha/pi)^3 factor
// from Remiddi-Laporta PLB 301 (1993) 440 Eq. (4)
double
amu_heavlept( const double mh )
{
  double t2, t4, t6, t8, t10, mh2, mh4, ln, zeta2;
  zeta2 = pi*pi/6;
  mh2 = mh*mh;
  mh4 = mh2*mh2;
  ln = log(mh);
  t2= (1.5*ZETA3-19.0/16)/mh2;
  t4= (-161.0/810*ln*ln - 16189.0/48600*ln + 13.0/18*ZETA3 - 161.0/9720*pi*pi - 831931.0/972000)/(mh2*mh2);
  t6= (17.0/36*ZETA3-13.0/224*zeta2-1840256147.0/3556224000-4381.0/120960*4.0*ln*ln-24761.0/317520*2.0*ln)/(mh4*mh2);
  t8= (7.0/20*ZETA3 - 2047.0/54000*zeta2 - 453410778211.0/1200225600000 - 5207.0/189000*4.0*ln*ln - 41940853.0/952560000*2.0*ln)/(mh4*mh4);
  t10=(5.0/18*ZETA3 - 1187.0/44550*zeta2 - 86251554753071.0/287550049248000 - 328337.0/14968800*4.0*ln*ln - 640572781.0/23051952000*2.0*ln)/(mh4*mh4*mh2);
  printf("# t2=%.11lg\n# t4=%.11lg\n# t6=%.11lg\n# t8=%.11lg\n# t10=%.11lg\n# sum= %.11lg\n", t2, t4, t6, t8, t10, t2+t4+t6+t8+t10);
  return(AlfQED*AlfQED*AlfQED/(pi*pi*pi)*(t2+t4+t6+t8+t10));
}

// contribution of a light lepton in a_mu, with ml=m_lepton/m_mu
// // from Remiddi-Laporta PLB 301 (1993) 440 Eq. (2)
double
amu_lightlept( const double ml )
{
  const double ml2 = ml*ml;
  const double ln = -log(ml);
  const double l2 = log(2.0);
  const double t1 = 2.0/3*pi*pi*ln;
  const double t2 = 59.0/270*pi*pi*pi*pi-3.0*ZETA3-10.0/3*pi*pi+2.0/3;
  const double t3 = ml*(4.0/3*pi*pi*ln-196.0/3*pi*pi*l2+424.0/9*pi*pi);
  const double t4 = ml2*(-2.0/3*ln*ln*ln+(pi*pi/9-20.0/3)*ln*ln-(16.0/135*pi*pi*pi*pi+4.0*ZETA3-32.0/9*pi*pi+61.0/3)*ln
	      +4.0/3*ZETA3*pi*pi-61.0/270*pi*pi*pi*pi+3.0*ZETA3+25.0/18*pi*pi-283.0/12);
  const double t5 = ml2*ml*(10.0/9*pi*pi*ln-11.0/9*pi*pi);
  const double t6 = ml2*ml2*(7.0/9*ln*ln*ln+41.0/18*ln*ln+13.0/9*pi*pi*ln+517.0/108*ln+ZETA3/2+191.0/216*pi*pi+13283.0/2592);
  printf("# t1=%.11lg\n# t2=%.11lg\n# t3=%.11lg\n# t4=%.11lg\n# t5=%.11lg\n# t6=%.11lg\n# sum= %.11lg\n", t1, t2, t3, t4, t5, t6, t1+t2+t3+t4+t5+t6);
  return(AlfQED*AlfQED*AlfQED/(pi*pi*pi)*(t1+t2+t3+t4+t5+t6));
}

// some kinda bessel function
void
bessk( const int n,
       const double x,
       double *bv )
{
  double bk, bkm, bkp, y, bi0, bi1, sq;
  const double tox=2.0/x;
  
  if(x<=2.0) {
    y=x/3.75;
    y*=y;
    bi0=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
	 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    bi1=x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
       +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));

    y=0.25*x*x;
    sq = log(x/2.0);
    bkm=(-sq*bi0)+(-0.57721566+y*(0.42278420
	 +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
	 +y*(0.10750e-3+y*0.74e-5))))));

    bk = (sq*bi1)+(1.0/x)*(1.0+y*(0.15443144
          +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
	  +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
    sq = (exp(-x)/sqrt(x));
    bkm=sq*(1.25331414+tox*(-0.7832358e-1
	 +tox*(0.2189568e-1+tox*(-0.1062446e-1+tox*(0.587872e-2
         +tox*(-0.251540e-2+tox*0.53208e-3))))));

    bk = sq*(1.25331414+tox*(0.23498619
	  +tox*(-0.3655620e-1+tox*(0.1504268e-1+tox*(-0.780353e-2
          +tox*(0.325614e-2+tox*(-0.68245e-3)))))));
  }

  bv[0]=bkm;
  bv[1]=bk;
  int j ;
  for(j=1;j<n;j++) {
    bkp=bkm+j*tox*bk;
    bkm=bk;
    bv[j+1]=bk=bkp;
  }
  return ;
}

// returns 1/4*tr(gamma_{idx[0]}*gamma_{idx[1]}...*gamma_{idx[n-1]})
int
diractrace( const int n ,
	    const int *idx)
{
  if(n&1)  return 0;
  if(n==0) return 1;
  if(n==2) return idx[0]==idx[1] ;

  int j, k, i, it ; 
  int idxred[n-2] ;
  for(it=0, i=0;i<n-1;i++) {
    if(idx[i]==idx[n-1]) {
      for(j=k=0;j<n-1;j++) {
	if(j!=i) idxred[k++] = idx[j];
      }
      // ternary-less version of i&1 ? -1 : 1 
      it += ( 1-2*(i&1) ) * diractrace(n-2,idxred);
    }
  }
  return it ;
}

// initialise the 8-index gamma combinations
void
init_g8( struct QED_kernel_temps *t )
{
  int i ;
  for( i = 0 ; i < 65536 ; i++ ) {
    // is in dta,rho,sig,mu,alf,nu,bet,lda order
    const int idx[8] = { (i/16384)&3 ,
			 (i/4096)&3 ,
			 (i/1024)&3 ,
			 (i/256)&3 ,
			 (i/64)&3 ,
			 (i/16)&3 ,
			 (i/4)&3 ,
			 i&3 } ;
    t -> G8[ i ] = 4*diractrace(8,idx) ;
  }
  return ;
}

void
ipihatFermLoop_antisym( const double x[4] ,
			const double y[4] ,
			const struct QED_kernel_temps t ,
			double vpihat[6][4][4][4] )
{
  const double mx[4] = { -x[0] , -x[1] , -x[2] , -x[3] } ;
  const double ymx[4] = { y[0]-x[0] , y[1]-x[1] , y[2]-x[2] , y[3]-x[3] } ;

  double vpihat1[4][4][4][4][4], vpihatr1[4][4][4][4];
  double vpihat2[4][4][4][4][4], vpihatr2[4][4][4][4];
  double vpihat3[4][4][4][4][4], vpihatr3[4][4][4][4];
  pihat1(x,y,t,vpihat1,vpihatr1);
  pihat1(ymx,mx,t,vpihat2,vpihatr2);
  pihat1(mx,ymx,t,vpihat3,vpihatr3);

  int rho, sig, mu, nu, lda ; 
  int rhosig=-1;
  for(rho=0;rho<4;rho++) {
    for(sig=rho+1;sig<4;sig++) {
      ++rhosig;
      for(mu=0;mu<4;mu++) {
	for(nu=0;nu<4;nu++) {
	  for(lda=0;lda<4;lda++) {
	    vpihat[rhosig][mu][nu][lda] =
	      0.5*( +vpihat1[rho][mu][nu][lda][sig]
		    +vpihat2[rho][nu][lda][mu][sig]
		    +x[rho]*vpihatr2[nu][lda][mu][sig]
		    +vpihat3[rho][lda][nu][mu][sig]
		    +x[rho]*vpihatr3[lda][nu][mu][sig]
		    - ( + vpihat1[sig][mu][nu][lda][rho]
			+ vpihat2[sig][nu][lda][mu][rho]
			+ x[sig]*vpihatr2[nu][lda][mu][rho]
			+ vpihat3[sig][lda][nu][mu][rho]
			+ x[sig]*vpihatr3[lda][nu][mu][rho] ));
	  }
	}
      }
    }
  }
  return ;
}

#define unroll_gamdta(gam,dta)			\
  Tmp = xab[alf][bta]*t.G8[ lexi_G8(alf,mu,bta,nu,gam,sig,dta,lda) ] ;	\
  t1[0] -= *fvP*Tmp ; fvP++ ;						\
  t1[1] -= *fvP*Tmp ; fvP++ ;						\
  t1[2] -= *fvP*Tmp ; fvP++ ;						\
  t1[3] -= *fvP*Tmp ; fvP++ ;						\
  tr1   -= *lP*Tmp ; lP++ ;						\

// sets pihat1 and pihatr1
void
pihat1( const double x[4] ,
	const double y[4] ,
	const struct QED_kernel_temps t ,
	double vpihat1[4][4][4][4][4],
	double vpihatr1[4][4][4][4] )
{
  const double pipi    = M_PI*M_PI ;
  const double twopipi = 2*pipi ;
  const double norm    = 8.2336243479292e-07; // 2/(2*pi)**8

  const double xmy[4] = { x[0]-y[0] , x[1]-y[1] , x[2]-y[2] , x[3]-y[3] } ;
  const double xn = v4norm(x);
  const double yn = v4norm(y);
  const double xmyn = v4norm(xmy);

  // bessel functions
  double bx[3], by[3], bxmy[3];
  bessk(2, xn, bx); bx[1] /= xn; bx[2] /= (xn*xn);
  bessk(2, yn, by); by[1] /= yn; by[2] /= (yn*yn);
  bessk(2, xmyn, bxmy); bxmy[1] /= xmyn; bxmy[2] /= (xmyn*xmyn);

  const double e[4] = { twopipi*y[0]*by[2] , twopipi*y[1]*by[2] ,
			twopipi*y[2]*by[2] , twopipi*y[3]*by[2] } ;
  const double f[4] = { pipi*y[0]*by[1] , pipi*y[1]*by[1] ,
			pipi*y[2]*by[1] , pipi*y[3]*by[1] } ;
  const double g[4] = { pipi*y[0]*by[0] , pipi*y[1]*by[0] ,
			pipi*y[2]*by[0] , pipi*y[3]*by[0] } ;

  // with h and ell written in these forms we save nearly two fp multiplies
  // per element of matrix h and matrix ell
  const double p = twopipi*by[0];
  const double q[4] = { 2*f[0] , 2*f[1] , 2*f[2] , 2*f[3] } ;
  // h is pi^2(y[rho].y[dta].b[1]-\delta_{rho,dta}b[0]).
  const double h[4][4] = { { f[0]*y[0] - 0.5*p , f[0]*y[1] , f[0]*y[2] , f[0]*y[3] } ,
			   { f[1]*y[0] , f[1]*y[1] - 0.5*p , f[1]*y[2] , f[1]*y[3] } ,
			   { f[2]*y[0] , f[2]*y[1] , f[2]*y[2] - 0.5*p , f[2]*y[3] } ,
			   { f[3]*y[0] , f[3]*y[1] , f[3]*y[2] , f[3]*y[3] - 0.5*p } } ;
  // ht is the transpose at it is more cache-coherent later on in the code
  const double ht[4][4] = { { h[0][0] , h[1][0] , h[2][0] , h[3][0] } ,
			    { h[0][1] , h[1][1] , h[2][1] , h[3][1] } ,
			    { h[0][2] , h[1][2] , h[2][2] , h[3][2] } ,
			    { h[0][3] , h[1][3] , h[2][3] , h[3][3] } } ;
  // ell is 2*\pi^2(y[rho].y[gam].b[2] - \delta_{rho,gam}b[1])
  const double s = twopipi*by[1];
  const double ell[4][4] = { { e[0]*y[0] - s , e[0]*y[1] , e[0]*y[2] , e[0]*y[3] } ,
			     { e[1]*y[0] , e[1]*y[1] - s , e[1]*y[2] , e[1]*y[3] } ,
			     { e[2]*y[0] , e[2]*y[1] , e[2]*y[2] - s , e[2]*y[3] } ,
			     { e[3]*y[0] , e[3]*y[1] , e[3]*y[2] , e[3]*y[3] - s } } ;

  // here we use the fact that fhatv is the same as h but with the
  // opposite sign of the delta so -delta + 2*delta = delta
  // and that 2*delta is in fact the variable "p"
  double fhatv[4][4], fv[4][4][4] ;
  int alf, bta, gam, dta, lda, sig, rho, mu, nu;
  for(rho=0;rho<4;rho++) {
    for(dta=0;dta<4;dta++) {
      fhatv[dta][rho] = h[rho][dta] ;
      for(gam=0;gam<4;gam++) {
	fv[gam][dta][rho] = pipi*(by[2]*y[gam]*y[dta]*y[rho]+
				  (+kro(rho,dta)*y[gam]
				   -kro(gam,rho)*y[dta]
				   -kro(gam,dta)*y[rho])*by[1]);	
      }
    }
    fhatv[rho][rho] += p ;
  }

  const double bx1bxmy1 = bx[1]*bxmy[1] ;
  const double bx1bxmy2 = bx[1]*bxmy[2] ;
  const double bx2bxmy1 = bx[2]*bxmy[1] ;
  const double bx2bxmy2 = bx[2]*bxmy[2] ;

  // precomputations
  double xab[4][4] , x2[4] , x3[4] , x4[4] , x5[4] ;
  for(alf=0;alf<4;alf++) {
    for(bta=0;bta<4;bta++) {
      xab[alf][bta] = x[alf]*xmy[bta]*bx2bxmy2 ;
    }
    x2[alf] = x[alf]*bx2bxmy1 ;
    x3[alf] = xmy[alf]*bx1bxmy2 ;
    x4[alf] = x[alf]*bx2bxmy1 ;
    x5[alf] = xmy[alf]*bx1bxmy2 ;
  }

  // flattened sigma,mu,nu,lda indexing
  int i ;
  for(i=0 ; i<256;i++ ) {
    sig = i/(64) ;
    mu  = (i/16)&3 ;
    nu  = (i/4)&3 ;
    lda = i&3 ;
      
    // initialize to zero our summations
    register double tr1 = 0.0 ;
    double t1[4] = { 0. , 0. , 0. , 0. } ;

    const double *fvP2 = (const double*)fv ;

    // alpha and beta are summation indices
    for(alf=0;alf<4;alf++) {

      const double *htP = (const double*)ht ;
      const double *fhatP = (const double*)fhatv ;
      
      for(bta=0;bta<4;bta++) {
	
	// this part of the code is hot so it gets unrolled
	register double Tmp ;
	const double *fvP = (const double*)fv, *lP = (const double*)ell ;
	unroll_gamdta(0,0);unroll_gamdta(0,1);unroll_gamdta(0,2);unroll_gamdta(0,3);
	unroll_gamdta(1,0);unroll_gamdta(1,1);unroll_gamdta(1,2);unroll_gamdta(1,3);
	unroll_gamdta(2,0);unroll_gamdta(2,1);unroll_gamdta(2,2);unroll_gamdta(2,3);
	unroll_gamdta(3,0);unroll_gamdta(3,1);unroll_gamdta(3,2);unroll_gamdta(3,3);
	      
	// the following was loop-fused as alf,bta,gam,dta are nothing
	// more than loop variables
	// again caching the multiplier in the loop that doesn't depend on rho
	const double tmp1 =
	  t.G8[ lexi_G6(alf,mu,bta,nu,sig,lda) ]*xab[alf][bta] ;
	const double tmp2 =
	  t.G8[ lexi_G6(alf,mu,nu,bta,sig,lda) ]*x2[alf] -
	  t.G8[ lexi_G6(mu,alf,nu,bta,sig,lda) ]*x3[alf] ;
	const double tmp4 =
	  t.G8[ lexi_G6(alf,mu,nu,sig,bta,lda) ]*x4[alf] -
	  t.G8[ lexi_G6(mu,alf,nu,sig,bta,lda) ]*x5[alf] ;
	const double tmp6 =
	  t.G8[ lexi_G6(mu,nu,alf,sig,bta,lda) ]*bx1bxmy1 ;

	t1[0] += -g[0]*tmp1 - *htP*(tmp2) - *fhatP*(tmp4) + *fvP2*tmp6 ; htP++ ; fhatP++ ; fvP2++ ;
	t1[1] += -g[1]*tmp1 - *htP*(tmp2) - *fhatP*(tmp4) + *fvP2*tmp6 ; htP++ ; fhatP++ ; fvP2++ ;
	t1[2] += -g[2]*tmp1 - *htP*(tmp2) - *fhatP*(tmp4) + *fvP2*tmp6 ; htP++ ; fhatP++ ; fvP2++ ;
	t1[3] += -g[3]*tmp1 - *htP*(tmp2) - *fhatP*(tmp4) + *fvP2*tmp6 ; htP++ ; fhatP++ ; fvP2++ ;
	
	tr1 -= p*tmp1 ;
	tr1 -= q[bta]*(tmp2+tmp4) ;
	tr1 += ell[alf][bta]*tmp6 ;
      }
    }
    // compute the last contribution and poke into vpihat1 and vpihatr1
    const double tmp = t.G8[ lexi_G4(mu,nu,sig,lda) ] * bx1bxmy1 ;
    for(rho=0;rho<4;rho++) {
      t1[rho] += g[rho]*tmp ;
      vpihat1[rho][mu][nu][lda][sig] = norm*t1[rho];
    }
    tr1 += p*tmp ;
    vpihatr1[mu][nu][lda][sig] = norm*tr1 ;
  }

  return ;
}

// vacuum polarisation
double
vac_pol( const double ml )
{
  const double a4 = 0.51747906167389938633;//PolyLog[4, 1/2]
  const double ml2 = ml*ml;
  const double ln = -log(ml);
  const double l2 = log(2.0);
  const double t1 = 2.0/9*ln*ln;
  const double t2 = (ZETA3-2.0/3*pi*pi*l2+pi*pi/9+31.0/27)*ln;
  const double t3 = 11.0/216*pi*pi*pi*pi-2.0/9*pi*pi*l2*l2-8.0/3*a4-1.0/9*l2*l2*l2*l2 -3.0*ZETA3+5.0/3*pi*pi*l2-25.0/18*pi*pi+1075.0/216;
  const double t4 = ml*(-13.0/18*pi*pi*pi-16.0/9*pi*pi*l2+3199.0/1080*pi*pi);
  const double t5 = ml2*(10.0/3*ln*ln-11.0/9*ln-14.0/3*pi*pi*l2-2.0*ZETA3+49.0/12*pi*pi-131.0/54);
  const double t6 = ml*ml2*(4.0/3*pi*pi*ln+35.0/12*pi*pi*pi-16.0/3*pi*pi*l2-5771.0/1080*pi*pi);
  return t1+t2+t3+t4+t5+t6 ; // should give (but it doesn't!) 20.9479242 for ml = 1.0/206.768262.
}

#ifdef pi
  #undef pi
#endif

#ifdef ZETA3
  #undef ZETA3
#endif
