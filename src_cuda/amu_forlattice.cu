/**
   @file amu_for_lattice.c
   @brief computes the QED kernel
*/
#include "KQED.h"      // definitions and what have you
#include "cu_util.h"

#include <vector>

// need to wrap constant vec arguments to CUDA kernels in a struct so they are
// passed by value (copied into constant global memory)
struct __attribute__((packed, aligned(8))) Vec4 {
  double x[4];
};
__host__
Vec4 vec4(const double xv[4]) {
  Vec4 pt;
  for (int i = 0; i < 4; ++i) {
    pt.x[i] = xv[i];
  }
  return pt;
}
struct __attribute__((packed, aligned(8))) OneKernel {
  double k[6][4][4][4] ;
};

#define CUDA_BLOCK_SIZE 32

// single point-like evaluation kernels
__global__
void
ker_QED_kernel_L0(
    const Vec4 *d_xv, const Vec4 *d_yv, unsigned n,
    const struct QED_kernel_temps t, OneKernel *d_kerv ) {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  QED_kernel_L0(d_xv[i].x, d_yv[i].x, t, (double (*)[4][4][4]) &d_kerv[i].k);
}
__global__
void
ker_pt_QED_kernel_L1(
    const Vec4 xv , const Vec4 yv ,
    const struct QED_kernel_temps t , double *d_kerv ) {
  QED_kernel_L1(xv.x, yv.x, t, (double (*)[4][4][4]) d_kerv);
}
__global__
void
ker_pt_QED_kernel_L2(
    const Vec4 xv , const Vec4 yv ,
    const struct QED_kernel_temps t , double *d_kerv ) {
  QED_kernel_L2(xv.x, yv.x, t, (double (*)[4][4][4]) d_kerv);
}
__global__
void
ker_pt_QED_Mkernel_L2(
    const double M, const Vec4 xv , const Vec4 yv ,
    const struct QED_kernel_temps t , double *d_kerv ) {
  QED_Mkernel_L2(M, xv.x, yv.x, t, (double (*)[4][4][4]) d_kerv);
}
__global__
void
ker_pt_QED_kernel_L3(
    const Vec4 xv , const Vec4 yv ,
    const struct QED_kernel_temps t , double *d_kerv ) {
  QED_kernel_L3(xv.x, yv.x, t, (double (*)[4][4][4]) d_kerv);
}
__global__
void
ker_pt_compute_all_kernels(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct QED_Kernels *K) {
  compute_all_kernels(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_all_kernels_SYMXY_v2(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct Kernels *K) {
  compute_all_kernels_SYMXY_v2(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_all_kernels_SYMXY(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct QED_Kernels *K) {
  compute_all_kernels_SYMXY(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_all_kernels_SYMXY0_v2(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct Kernels *K) {
  compute_all_kernels_SYMXY0_v2(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_all_kernels_SYMXY0(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct QED_Kernels *K) {
  compute_all_kernels_SYMXY0(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_all_Mkernels(
    const Vec4 M, const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct QED_Kernels *K) {
  compute_all_Mkernels(M.x, xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_all_Mkernels_v2(
    const Vec4 M, const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct QED_Kernels *K) {
  compute_all_Mkernels_v2(M.x, xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_con_kernels(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct QED_Kernels *K) {
  compute_con_kernels(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_con_kernels_v2(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct Kernels *K) {
  compute_con_kernels_v2(xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_con_kernelsM_L2(
    const Vec4 M, const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct Kernels *K) {
  compute_con_kernelsM_L2(M.x, xv.x, yv.x, t, K);
}
__global__
void
ker_pt_compute_sub_kernelsM_L2(
    const Vec4 M, const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, struct Kernels *K) {
  compute_sub_kernelsM_L2(M.x, xv.x, yv.x, t, K);
}
__global__
void
ker_pt_ipihatFermLoop_antisym(
    const Vec4 xv, const Vec4 yv,
    const struct QED_kernel_temps t, double *d_vpihat) {
  ipihatFermLoop_antisym(xv.x, yv.x, t, (double (*)[4][4][4]) d_vpihat);
}
__global__
void
ker_ipihatFermLoop_antisym(
    const Vec4 *d_xv, const Vec4 *d_yv, unsigned n,
    const struct QED_kernel_temps t, OneKernel *d_vpihat) {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  ipihatFermLoop_antisym(d_xv[i].x, d_yv[i].x, t, (double (*)[4][4][4]) &d_vpihat[i].k);
}


__host__
void
cu_pt_QED_kernel_L0(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , double kerv[6][4][4][4] ) {
  // double *d_kerv;
  OneKernel *d_kerv;
  size_t sizeof_kerv = 6*4*4*4*sizeof(double);
  checkCudaErrors(cudaMalloc(&d_kerv, sizeof_kerv));
  Vec4 *d_xv, *d_yv;
  checkCudaErrors(cudaMalloc(&d_xv, sizeof(Vec4)));
  checkCudaErrors(cudaMalloc(&d_yv, sizeof(Vec4)));
  checkCudaErrors(cudaMemcpy(d_xv, xv, sizeof(Vec4), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_yv, yv, sizeof(Vec4), cudaMemcpyHostToDevice));
  ker_QED_kernel_L0<<<1,1>>>( d_xv, d_yv, 1, t, d_kerv );
  checkCudaErrors(cudaMemcpy(kerv, d_kerv, sizeof_kerv, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_kerv));
  checkCudaErrors(cudaFree(d_xv));
  checkCudaErrors(cudaFree(d_yv));
}
__host__
void
cu_arr_QED_kernel_L0(
    const Vec4 *xv, const Vec4 *yv, const unsigned n,
    const struct QED_kernel_temps t , OneKernel *kerv ) {
  Vec4 *d_xv, *d_yv;
  checkCudaErrors(cudaMalloc(&d_xv, n*sizeof(Vec4)));
  checkCudaErrors(cudaMalloc(&d_yv, n*sizeof(Vec4)));
  checkCudaErrors(cudaMemcpy(d_xv, xv, n*sizeof(Vec4), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_yv, yv, n*sizeof(Vec4), cudaMemcpyHostToDevice));
  OneKernel *d_kerv;
  size_t sizeof_kerv = n*sizeof(OneKernel);
  checkCudaErrors(cudaMalloc(&d_kerv, sizeof_kerv));
  unsigned n_threads = CUDA_BLOCK_SIZE;
  unsigned n_blocks = (n + n_threads - 1) / n_threads;
  ker_QED_kernel_L0<<<n_blocks,n_threads>>>(d_xv, d_yv, n, t, d_kerv);
  checkCudaErrors(cudaMemcpy(kerv, d_kerv, sizeof_kerv, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_xv));
  checkCudaErrors(cudaFree(d_yv));
  checkCudaErrors(cudaFree(d_kerv));
}
__host__
void
cu_pt_QED_kernel_L1(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , double kerv[6][4][4][4] ) {
  double *d_kerv;
  size_t sizeof_kerv = 6*4*4*4*sizeof(double);
  checkCudaErrors(cudaMalloc(&d_kerv, sizeof_kerv));
  ker_pt_QED_kernel_L1<<<1,1>>>( vec4(xv), vec4(yv), t, d_kerv );
  checkCudaErrors(cudaMemcpy(kerv, d_kerv, sizeof_kerv, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_kerv));
}
__host__
void
cu_pt_QED_kernel_L2(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , double kerv[6][4][4][4] ) {
  double *d_kerv;
  size_t sizeof_kerv = 6*4*4*4*sizeof(double);
  checkCudaErrors(cudaMalloc(&d_kerv, sizeof_kerv));
  ker_pt_QED_kernel_L2<<<1,1>>>( vec4(xv), vec4(yv), t, d_kerv );
  checkCudaErrors(cudaMemcpy(kerv, d_kerv, sizeof_kerv, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_kerv));
}
__host__
void
cu_pt_QED_Mkernel_L2(
    const double M, const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , double kerv[6][4][4][4] ) {
  double *d_kerv;
  size_t sizeof_kerv = 6*4*4*4*sizeof(double);
  checkCudaErrors(cudaMalloc(&d_kerv, sizeof_kerv));
  ker_pt_QED_Mkernel_L2<<<1,1>>>( M, vec4(xv), vec4(yv), t, d_kerv );
  checkCudaErrors(cudaMemcpy(kerv, d_kerv, sizeof_kerv, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_kerv));
}
__host__
void
cu_pt_QED_kernel_L3(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , double kerv[6][4][4][4] ) {
  double *d_kerv;
  size_t sizeof_kerv = 6*4*4*4*sizeof(double);
  checkCudaErrors(cudaMalloc(&d_kerv, sizeof_kerv));
  ker_pt_QED_kernel_L3<<<1,1>>>( vec4(xv), vec4(yv), t, d_kerv );
  checkCudaErrors(cudaMemcpy(kerv, d_kerv, sizeof_kerv, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_kerv));
}

__host__
void
cu_pt_compute_all_kernels(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct QED_Kernels *K ) {
  struct QED_Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(QED_Kernels)));
  ker_pt_compute_all_kernels<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(QED_Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_all_kernels_SYMXY_v2(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct Kernels *K ) {
  struct Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(Kernels)));
  ker_pt_compute_all_kernels_SYMXY_v2<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_all_kernels_SYMXY(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct QED_Kernels *K ) {
  struct QED_Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(QED_Kernels)));
  ker_pt_compute_all_kernels_SYMXY<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(QED_Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_all_kernels_SYMXY0_v2(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct Kernels *K ) {
  struct Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(Kernels)));
  ker_pt_compute_all_kernels_SYMXY0_v2<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_all_kernels_SYMXY0(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct QED_Kernels *K ) {
  struct QED_Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(QED_Kernels)));
  ker_pt_compute_all_kernels_SYMXY0<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(QED_Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_all_Mkernels(
    const double M[4], const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct QED_Kernels *K ) {
  struct QED_Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(QED_Kernels)));
  ker_pt_compute_all_Mkernels<<<1,1>>>( vec4(M), vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(QED_Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_all_Mkernels_v2(
    const double M[4], const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct QED_Kernels *K ) {
  struct QED_Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(QED_Kernels)));
  ker_pt_compute_all_Mkernels_v2<<<1,1>>>( vec4(M), vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(QED_Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_con_kernels(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct QED_Kernels *K ) {
  struct QED_Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, sizeof(QED_Kernels)));
  ker_pt_compute_con_kernels<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, sizeof(QED_Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_con_kernels_v2(
    const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct Kernels K[3] ) {
  struct Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, 3*sizeof(Kernels)));
  ker_pt_compute_con_kernels_v2<<<1,1>>>( vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, 3*sizeof(Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_con_kernelsM_L2(
    const double M[4], const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct Kernels K[3] ) {
  struct Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, 3*sizeof(Kernels)));
  ker_pt_compute_con_kernelsM_L2<<<1,1>>>( vec4(M), vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, 3*sizeof(Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_pt_compute_sub_kernelsM_L2(
    const double M[4], const double xv[4] , const double yv[4] ,
    const struct QED_kernel_temps t , struct Kernels K[3] ) {
  struct Kernels *d_K;
  checkCudaErrors(cudaMalloc(&d_K, 3*sizeof(Kernels)));
  ker_pt_compute_sub_kernelsM_L2<<<1,1>>>( vec4(M), vec4(xv), vec4(yv), t, d_K );
  checkCudaErrors(cudaMemcpy(K, d_K, 3*sizeof(Kernels), cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_K));
}
__host__
void
cu_arr_ipihatFermLoop_antisym(
    const Vec4 *xv , const Vec4 *yv , const unsigned n ,
    const struct QED_kernel_temps t , OneKernel *vpihat ) {
  Vec4 *d_xv, *d_yv;
  checkCudaErrors(cudaMalloc(&d_xv, n*sizeof(Vec4)));
  checkCudaErrors(cudaMalloc(&d_yv, n*sizeof(Vec4)));
  checkCudaErrors(cudaMemcpy(d_xv, xv, n*sizeof(Vec4), cudaMemcpyHostToDevice));
  checkCudaErrors(cudaMemcpy(d_yv, yv, n*sizeof(Vec4), cudaMemcpyHostToDevice));
  OneKernel *d_vpihat;
  size_t sizeof_vpihat = n*sizeof(OneKernel);
  checkCudaErrors(cudaMalloc(&d_vpihat, sizeof_vpihat));
  unsigned n_threads = CUDA_BLOCK_SIZE;
  unsigned n_blocks = (n + n_threads - 1) / n_threads;
  ker_ipihatFermLoop_antisym<<<n_blocks,n_threads>>>( d_xv, d_yv, n, t, d_vpihat );
  checkCudaErrors(cudaMemcpy(vpihat, d_vpihat, sizeof_vpihat, cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaFree(d_xv));
  checkCudaErrors(cudaFree(d_yv));
  checkCudaErrors(cudaFree(d_vpihat));
}


// pretend discrete volume
static const size_t LVOLUME = (2*2*2*4) ;

__host__
static void
HLBL_crds( double x[4] ,
	   const size_t idx )
{      
  x[0] = idx%2      ; x[0] = ( x[0] > 1 ) ? x[0] - 2 : x[0] ;
  x[1] = (idx/2)%2  ; x[1] = ( x[1] > 1 ) ? x[1] - 2 : x[1] ;
  x[2] = (idx/4)%2  ; x[2] = ( x[2] > 1 ) ? x[2] - 2 : x[2] ;
  x[3] = (idx/8) ; x[3] = ( x[3] > 2 ) ? x[3] - 4 : x[3] ;
  return ;
}

// test that all elements of kernel are zero
__host__
static bool
kernel_nonzero( const double *kerv )
{
  int i ;
  for( i = 0 ; i < 384 ; i++ ) {
    if( fabs( *kerv ) > 1E-14 ) {
      return true ;
    }
    kerv++ ;
  }
  return false ;
}

__host__
static void
printbar( void )
{
  fprintf( stdout , "\n-----------------------------------------------------\n" ) ;
}

__host__
static void
print_header( int i )
{
  printbar() ;
  fprintf( stdout , "=== EXAMPLE %d ===\n", i ) ;
}

// printing utility
__host__
static void
print_kernels( const double *kerv ,
	       const double *kerv2 ,
	       FILE *file )
{
  printbar() ;
  
  int i ;
  for( i = 0 ; i < 384 ; i++ ) {
    const int rhosig = i/(64) ;
    const int mu = (i/16)&3 ;
    const int nu = (i/4)&3 ;
    const int lambda = i&3 ;

    fprintf( file , "rhosig= %d\tmu= %d\tnu= %d\tlda= %d\t%lg   %lg\n",
	     rhosig, mu, nu, lambda, (*kerv) , (*kerv2 ) ) ;
    kerv++ ; kerv2++ ;
  }
  printbar() ;
  return ;
}

// contracts the kernel with pihat
__host__
static void
example1( const struct QED_kernel_temps t )
{  
  static const double Mv = 2.0 ; // ratio of lepton-loop mass to muon mass
  const double pref = (pow(4.0*M_PI*AlfQED,3)/3)*(2.0*M_PI*M_PI)*(4*M_PI);
  const double convf = pow(Mv,7);
  const double norm_fermloop = pow(1.0,4);
  const double fnorm = 2.0*pref*convf*norm_fermloop ;

  start_timer() ;

  std::vector<Vec4> xvs, yvs, xvMvs, yvMvs;
  
  double x,y,z ;
  for( x = 0.005 ; x < 1 ; x += 0.2 ) {
    for( y = 0.005 ; y < M_PI ; y += M_PI/4. ) {
      for( z = 0.005 ; z < 1 ; z += 0.2 ) { 

	const double vi[3] = { x , y , z } ;

	const double co = cos(vi[1]);
	const double si = sin(vi[1]);
    
	const double xv[4] = { 0 , 0 , vi[0]*si , vi[0]*co } ;
    	const double yv[4] = { 0 , 0 , 0 , vi[2] } ;

	const double yvMv[4] = { 0 , 0 , 0 , vi[2]*Mv } ; 
	const double xvMv[4] = { 0 , 0 , vi[0]*si*Mv , vi[0]*co*Mv } ;

        xvs.push_back(vec4(xv));
        yvs.push_back(vec4(yv));
        xvMvs.push_back(vec4(xvMv));
        yvMvs.push_back(vec4(yvMv));
      }
    }
  }

  std::vector<OneKernel> kervs(xvs.size());
  std::vector<OneKernel> pihats(xvs.size());
  cu_arr_QED_kernel_L0(
      xvs.data(), yvs.data(), xvs.size(), t, kervs.data());
  cu_arr_ipihatFermLoop_antisym(
      xvMvs.data(), yvMvs.data(), xvs.size(), t, pihats.data());

  unsigned i = 0;
  for( x = 0.005 ; x < 1 ; x += 0.2 ) {
    for( y = 0.005 ; y < M_PI ; y += M_PI/4. ) {
      for( z = 0.005 ; z < 1 ; z += 0.2 ) { 

	const double vi[3] = { x , y , z } ;
	const double si = sin(vi[1]);

	const double *pi = (const double*)&pihats[i].k ;
	const double *kp = (const double*)&kervs[i].k ;
	double tmp = 0.0;
	int idx ;
	for( idx = 0 ; idx < 384 ; idx++ ) {
	  tmp += *pi * ( *kp ) ;
	  pi++ ; kp++ ;
	}
	tmp *= fnorm;
	fprintf( stdout , "x= %lf beta= %lf  y= %lf  iPihat*L_QED= %.11lg\n",
		 vi[0], vi[1], vi[2], 2*pow(vi[0],3.0)*pow(vi[2],4.0)*si*si*tmp);

        i++;
      }
    }
  }
  print_time() ;

  return ;
}

// computes the L0 kernel with xy and with x = 0
__host__
static void
example2( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example2.txt" , "w" ) ;
  
  double kerv1[6][4][4][4] KQED_ALIGN ; 
  double kerv2[6][4][4][4] KQED_ALIGN ;

  const double zero[4] = { 0. , 0. , 0. , 0. } ;
  const double vi[3] = { 1.228 , 0.209243 , 0.0484223 } ;

  const double co = cos(vi[1]);
  const double si = sin(vi[1]);
  
  const double yv[4] = { 0 , 0 , 0 , vi[2] } ; 
  const double xv[4] = { 0 , 0 , vi[0]*si , vi[0]*co } ;

  start_timer() ;
  
  printf("# vi[0]= %lg   vi[1] = %lg   vi[2] = %lg\n",
	 vi[0] , vi[1] , vi[2] ) ;

  cu_pt_QED_kernel_L0( xv , yv , t , kerv1 ) ;
  cu_pt_QED_kernel_L0( zero , yv , t , kerv2 ) ;
  
  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  cu_pt_QED_kernel_L0( yv , xv , t , kerv1 ) ;
  cu_pt_QED_kernel_L0( xv , zero , t , kerv2 ) ;
  
  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  print_time() ;

  fclose( file ) ;
  
  return ;
}

// tests some identities for the different L-kernels
__host__
void
example3( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example3.txt" , "w" ) ;
  double kerv1[6][4][4][4] KQED_ALIGN , kerv2[6][4][4][4] KQED_ALIGN ;
  double kerv3[6][4][4][4] KQED_ALIGN , kerv4[6][4][4][4] KQED_ALIGN ;
  const double zero[4] = { 0 , 0 , 0 , 0 } ; 
  const double xv[4] = { 1 , 2 , 3 , 4 } ;
  const double yv[4] = { 4 , 3 , 2 , 1 } ;
  start_timer() ;
  printbar() ;
  // test that L_0(0,0) == 0 and that L_1(x,x) == 0
  fprintf( stdout , "\nTesting that L_0(0,0) = L_1(x,x) == 0\n" ) ;
  cu_pt_QED_kernel_L0( zero , zero , t , kerv1 ) ;
  cu_pt_QED_kernel_L1( xv , xv , t , kerv2 ) ;
  if( kernel_nonzero( (const double*)kerv1 ) ||
      kernel_nonzero( (const double*)kerv2 ) ) {
    print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;
    fprintf( stderr , "!Test Failed!\n" ) ;
  }
  // test that L_2(0,y) == 0 and that L_2(x,0) == 0
  fprintf( stdout , "\nTesting that L_2(0,y) = L_2(x,0) == 0\n" ) ;
  cu_pt_QED_kernel_L2( zero , yv , t , kerv1 ) ;
  cu_pt_QED_kernel_L2( xv , zero , t , kerv2 ) ;
  if( kernel_nonzero( (const double*)kerv1 ) ||
      kernel_nonzero( (const double*)kerv2 ) ) {
    print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;
    fprintf( stderr , "\n!Test Failed!\n" ) ;
  }
  // test that L_3(x,x) = L_3(0,y) == L_3(y,y) == L_3(0,x) == 0
  fprintf( stdout , "\nTesting that L_3(x,x) = L_3(0,y) == L_3(y,y) == L_3(0,x) ==  0\n" ) ;
  cu_pt_QED_kernel_L3( xv , xv , t , kerv1 ) ;
  cu_pt_QED_kernel_L3( zero , yv , t , kerv2 ) ;
  cu_pt_QED_kernel_L3( zero , xv , t , kerv3 ) ;
  cu_pt_QED_kernel_L3( yv , yv , t , kerv4 ) ;
  if( kernel_nonzero( (const double*)kerv1 ) ||
      kernel_nonzero( (const double*)kerv2 ) ||
      kernel_nonzero( (const double*)kerv3 ) ||
      kernel_nonzero( (const double*)kerv4 ) ) {
    print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;
    print_kernels( (const double*)kerv3 , (const double*)kerv4 , file ) ;
    fprintf( stderr , "\n!Test Failed!\n" ) ;
  }
  printbar() ;
  fclose( file ) ;
  print_time() ;
  return ;
}

// computes the L0 kernel triggering the MYSWAP
__host__
static void
example4( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example4.txt" , "w" ) ;
  
  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;

  const double x[4] = { 1 , 1 , 1 , 1.1 } ;
  const double y[4] = { 1 , 1 , 1 , 1 } ;

  start_timer() ;

  cu_pt_QED_kernel_L0( x , y , t , kerv1 ) ;

  cu_pt_QED_kernel_L0( y , x , t , kerv2 ) ; 

  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  print_time() ;

  fclose( file ) ;
  return ;
}

// computes the L0 kernel triggering the taylor expansions
__host__
static void
example5( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example5.txt" , "w" ) ;
  
  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;

  const double x[4] = { 0.001 , 0.001 , 0.001 , 0.001 } ;
  const double y[4] = { 1 , 1 , 1 , 1 } ;

  start_timer() ;

  cu_pt_QED_kernel_L0( x , y , t , kerv1 ) ;

  cu_pt_QED_kernel_L0( y , x , t , kerv2 ) ; 

  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  print_time() ;

  fclose( file ) ;
  
  return ;
}

// computes the L0 kernel for some realistic example
__host__
static void
example6( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example6.txt" , "w" ) ;
  
  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;

  const double fac = 0.046241301098278 ; 
  const double x[4] = { 1*fac , 3*fac , 3*fac , 52*fac };
  const double y[4] = { 1*fac , 1*fac , 1*fac, 1*fac };
  
  start_timer() ;

  cu_pt_QED_kernel_L0( x , y , t , kerv1 ) ;

  cu_pt_QED_kernel_L0( y , x , t , kerv2 ) ; 

  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  print_time() ;

  fclose( file ) ;
  
  return ;
}

// computes the L0 kernel triggering taylorx and xmyswap routines
static void
example7( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example7.txt" , "w" ) ;
  
  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;

  const double x[4] = { 0.001 , 0.001 , 0.001 , 0.001 };
  const double y[4] = { 0.001 , 0.001 , 0.001 , -0.001 };
  
  start_timer() ;

  cu_pt_QED_kernel_L0( x , y , t , kerv1 ) ;

  cu_pt_QED_kernel_L0( y , x , t , kerv2 ) ; 

  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  print_time() ;

  fclose( file ) ;
  
  return ;
}

// see if we can use L(x,0) instead of L(x,x)
__host__
static void
example9( const struct QED_kernel_temps t )
{
  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;

  start_timer() ;

  printbar() ;
  fprintf( stdout , "Testing that L(x,x) = -L_{mu<->lambda}(x,0)\n" ) ;
  
  double theta ;
  for( theta = 0.1 ; theta < 1 ; theta+= 0.05 ) {
  
    const double zero[4] = { 0,0,0,0 };
    const double x[4] = { theta,2*theta,3*theta,4*theta };
    
    cu_pt_QED_kernel_L0( x , x , t , kerv1 ) ;
    cu_pt_QED_kernel_L0( x , zero , t , kerv2 ) ;
    
    size_t i , j , k , l ;
    for( i = 0 ; i < 6 ; i++ ) {
      for( j = 0 ; j < 4 ; j++ ) {
	for( k = 0 ; k < 4 ; k++ ) {
	  for( l = 0 ; l < 4 ; l++ ) {
	    if( fabs( kerv1[i][j][k][l] + kerv2[i][l][k][j] ) > 1E-14 ) {
	      fprintf( stderr , "Problem (%zu %zu %zu %zu) %e != %e\n" ,
		       i , j , k , l ,
		       kerv1[i][j][k][l] , -kerv2[i][l][k][j] ) ;
	    }
	  }
	}
      }
    }
  }
  print_time() ;
  printbar() ;

  return ;
}

// check L(x,0) and L(0,x)
__host__
static void
example10( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example10.txt" , "w" ) ;

  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;

  const double zero[4] = { 0., 0., 0., 0. };
  const double x[4] = { 0.1 ,0.1 ,0.1 ,0.1 };
  
  start_timer() ;

  cu_pt_QED_kernel_L0( zero , x , t , kerv1 ) ;
  cu_pt_QED_kernel_L0( x , zero , t , kerv2 ) ;

  print_kernels( (const double*)kerv1 , (const double*)kerv2 , file ) ;

  print_time() ;
  fclose( file ) ;

  return ;
}

// test that L(x,y) == -L(-x,-y)
__host__
static void
example11( const struct QED_kernel_temps t )
{
  FILE *file = fopen( "example11.txt" , "w" ) ;

  double kerv1[6][4][4][4] KQED_ALIGN ;
  double kerv2[6][4][4][4] KQED_ALIGN ;
  
  start_timer() ;

  size_t i , j ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    double y[4] ;
    HLBL_crds( y , i ) ;
    const double my[4] = {-y[0],-y[1],-y[2],-y[3]} ;

    for( j = 0 ; j < LVOLUME ; j++ ) {

      double x[4] ;
      HLBL_crds( x , j ) ;
      const double mx[4] = {-x[0],-x[1],-x[2],-x[3]} ;

      
      size_t n ;
      for( n = 0 ; n < 4 ; n++ ) {
	
	switch( n ) {
	case 0 :
	  cu_pt_QED_kernel_L0( x , y , t , kerv1 ) ;
	  cu_pt_QED_kernel_L0( mx , my , t , kerv2 ) ;
	  break ;
	case 1 :
	  cu_pt_QED_kernel_L1( x , y , t , kerv1 ) ;
	  cu_pt_QED_kernel_L1( mx , my , t , kerv2 ) ;
	  break ;
	case 2 :
	  cu_pt_QED_kernel_L2( x , y , t , kerv1 ) ;
	  cu_pt_QED_kernel_L2( mx , my , t , kerv2 ) ;
	  break ;
	case 3 :
	  cu_pt_QED_kernel_L3( x , y , t , kerv1 ) ;
	  cu_pt_QED_kernel_L3( mx , my , t , kerv2 ) ;
	  break ;
	default :
	  break ;
	}
	
	const double *p1 = (const double*)kerv1 ;
	const double *p2 = (const double*)kerv2 ;
	
	size_t i ;
	bool failure = false ;
	for( i = 0 ; i < 384 ; i++ ) {
	  if( fabs( *p1 + *p2 ) > 1E-14 ) {
	    fprintf( stderr , "Problem L_%zu (%zu) %e != %e\n" ,
		     n , i , *p1 , -*p2 ) ;
	    p1++ ; p2++ ;
	    failure = true ;
	  }
	}

	if( failure ) {
	  fprintf( stdout , "Test L^(%zu)(x,y) = -L^(%zu)(-x,-y) failed\n" ,
		   n , n ) ;
	}
	
      }// loop on kernels
    }
  }
  
  print_time() ;
  printbar() ;
  fclose( file ) ;

  return ;
}

__host__
static bool
equivalent_kernels( double K1[6][4][4][4] ,
		    double K2[6][4][4][4] )
{
  const double *pK1 = (const double*)K1 ;
  const double *pK2 = (const double*)K2 ;
  size_t i ;
  for( i = 0 ; i < 384 ; i++ ) {
    if( fabs( *pK1 - *pK2 ) > 1E-14 ) {
      return false ;
    }
    pK1++ ; pK2++ ;
  }
  return true ;
}

// test that our compute_all_kernels routine works
__host__
static void
example12( const struct QED_kernel_temps t )
{
  struct QED_Kernels K ;

  printbar() ;
  fprintf( stdout , "Checking compute_all_kernels routine\n" ) ;

  start_timer() ;
  
  size_t i , j ;
  double x[4] , y[4] ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    HLBL_crds( y , i ) ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      HLBL_crds( x , j ) ;
  
      cu_pt_compute_all_kernels( x , y , t , &K ) ;
      
      double Ltest[6][4][4][4] KQED_ALIGN ;
      size_t tests_failed = 0 ;
  
      // check L1 kernels
      cu_pt_QED_kernel_L1( x , y , t , Ltest ) ;
      if( !equivalent_kernels( Ltest , K.L1.xy ) ) {
	fprintf( stderr , "Inequivalent kernels L1xy\n" ) ;
	tests_failed++ ;
      }
      cu_pt_QED_kernel_L1( y , x , t , Ltest ) ;
      if( !equivalent_kernels( Ltest , K.L1.yx ) ) {
	fprintf( stderr , "Inequivalent kernels L1yx\n" ) ;
	tests_failed++ ;
      }
      // check L2 kernels
      cu_pt_QED_kernel_L2( x , y , t , Ltest ) ;
      if( !equivalent_kernels( Ltest , K.L2.xy ) ) {
	fprintf( stderr , "Inequivalent kernels L2xy\n" ) ;
	tests_failed++ ;
      }
      cu_pt_QED_kernel_L2( y , x , t , Ltest ) ;
      if( !equivalent_kernels( Ltest , K.L2.yx ) ) {
	fprintf( stderr , "Inequivalent kernels L2yx\n" ) ;
	tests_failed++ ;
      }
      // check L3 kernels
      cu_pt_QED_kernel_L3( x , y , t , Ltest ) ;
      if( !equivalent_kernels( Ltest , K.L3.xy ) ) {
	fprintf( stderr , "Inequivalent kernels L3xy\n" ) ;
	tests_failed++ ;
      }
      cu_pt_QED_kernel_L3( y , x , t , Ltest ) ;
      if( !equivalent_kernels( Ltest , K.L3.yx ) ) {
	fprintf( stderr , "Inequivalent kernels L3yx\n" ) ;
	tests_failed++ ;
      }
      
      if( tests_failed != 0 ) {
	fprintf( stderr , "%zu tests failed\n" , tests_failed ) ;
	fprintf( stderr , "x (%f,%f,%f,%f) y (%f,%f,%f,%f)\n" ,
		 x[0],x[1],x[2],x[3],y[0],y[1],y[2],y[3] ) ;
      }
    }
  }
  printbar() ;
  print_time() ;
  
  return ;
}

// test our SYMXY compute all kernels code
__host__
static void
example13( const struct QED_kernel_temps t )
{
  struct Kernels K1 , K2 ;

  printbar() ;
  fprintf( stdout , "Checking compute_all_kernels_SYMXY routine\n" ) ;

  start_timer() ;
  size_t X , Y ;
  double x[4] , y[4] ;
  for( X = 0 ; X < LVOLUME ; X++ ) {
    HLBL_crds( y , X ) ;
    for( Y = 0 ; Y < LVOLUME ; Y++ ) {
      HLBL_crds( x , Y ) ;
  
      cu_pt_compute_all_kernels_SYMXY_v2( x , y , t , &K1 ) ;
      cu_pt_compute_all_kernels_SYMXY_v2( y , x , t , &K2 ) ;

      size_t i , j , k , l ;
      for( i = 0 ; i < 6 ; i++ ) {
	for( j = 0 ; j < 4 ; j++ ) {
	  for( k = 0 ; k < 4 ; k++ ) {
	    for( l = 0 ; l < 4 ; l++ ) {
	      // L0
	      if( fabs( K1.L0[i][j][k][l] - K2.L0[i][k][j][l] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY L0 failure (%zu %zu %zu %zu) %e != %e\n" ,
			 i , j , k , l , K1.L0[i][j][k][l] , K2.L0[i][k][j][l] ) ; 
	      }
	      // L1
	      if( fabs( K1.L1[i][j][k][l] - K2.L1[i][k][j][l] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY L1 failure (%zu %zu %zu %zu) %e != %e\n" ,
			 i , j , k , l , K1.L1[i][j][k][l] , K2.L1[i][k][j][l] ) ; 
	      }
	      // L2
	      if( fabs( K1.L2[i][j][k][l] - K2.L2[i][k][j][l] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY L2 failure (%zu %zu %zu %zu) %e != %e\n" ,
			 i , j , k , l , K1.L2[i][j][k][l] , K2.L2[i][k][j][l] ) ; 
	      }
	      // L3
	      if( fabs( K1.L3[i][j][k][l] - K2.L3[i][k][j][l] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY L3 failure (%zu %zu %zu %zu) %e != %e\n" ,
			 i , j , k , l , K1.L3[i][j][k][l] , K2.L3[i][k][j][l] ) ; 
	      }
	      //
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  return ;
}

// test our SYMXY0 compute all kernels code
__host__
static void
example14( const struct QED_kernel_temps t )
{
  struct Kernels K1xy , K1yx , K2xy , K2yx , K3xy , K3yx ;

  printbar() ;
  fprintf( stdout , "Checking compute_all_kernels_SYMXY0 routine\n" ) ;

  start_timer() ;
  size_t X , Y ;
  double x[4] , y[4] ;
  for( X = 0 ; X < LVOLUME ; X++ ) {
    HLBL_crds( y , X ) ;
    for( Y = 0 ; Y < LVOLUME ; Y++ ) {
      HLBL_crds( x , Y ) ;
  
      cu_pt_compute_all_kernels_SYMXY0_v2( x , y , t , &K1xy ) ;
      cu_pt_compute_all_kernels_SYMXY0_v2( y , x , t , &K1yx ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;
      const double ymx[4] = { -xmy[0] , -xmy[1] , -xmy[2] , -xmy[3] } ;
      cu_pt_compute_all_kernels_SYMXY0_v2( x , xmy , t , &K2xy ) ;
      cu_pt_compute_all_kernels_SYMXY0_v2( xmy , x , t , &K2yx ) ;
      cu_pt_compute_all_kernels_SYMXY0_v2( y , ymx , t , &K3xy ) ;
      cu_pt_compute_all_kernels_SYMXY0_v2( ymx , y , t , &K3yx ) ;

      size_t i , j , k , l ;
      for( i = 0 ; i < 6 ; i++ ) {
	for( j = 0 ; j < 4 ; j++ ) {
	  for( k = 0 ; k < 4 ; k++ ) {
	    for( l = 0 ; l < 4 ; l++ ) {
	      // L0
	      if( fabs( K1xy.L0[i][j][k][l] - K1yx.L0[i][k][j][l] ) > 1E-12 ||
		  fabs( K1xy.L0[i][j][k][l] + K2xy.L0[i][l][k][j] ) > 1E-12 ||
		  fabs( K1xy.L0[i][j][k][l] + K2yx.L0[i][k][l][j] ) > 1E-12 ||
		  fabs( K1xy.L0[i][j][k][l] + K3xy.L0[i][l][j][k] ) > 1E-12 ||
		  fabs( K1xy.L0[i][j][k][l] + K3yx.L0[i][j][l][k] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY0 L0 failure (%zu %zu %zu %zu) %e != %e "
			 "!= %e != %e != %e != %e\n" ,  
			 i , j , k , l ,
			 +K1xy.L0[i][j][k][l] , +K1yx.L0[i][k][j][l] ,
			 -K2xy.L0[i][l][k][j] , -K2yx.L0[i][k][l][j] ,
			 -K3xy.L0[i][l][j][k] , -K3yx.L0[i][j][l][k] ) ; 
	      }
	      // L1
	      if( fabs( K1xy.L1[i][j][k][l] - K1yx.L1[i][k][j][l] ) > 1E-12 ||
		  fabs( K1xy.L1[i][j][k][l] + K2xy.L1[i][l][k][j] ) > 1E-12 ||
		  fabs( K1xy.L1[i][j][k][l] + K2yx.L1[i][k][l][j] ) > 1E-12 ||
		  fabs( K1xy.L1[i][j][k][l] + K3xy.L1[i][l][j][k] ) > 1E-12 ||
		  fabs( K1xy.L1[i][j][k][l] + K3yx.L1[i][j][l][k] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY0 L1 failure (%zu %zu %zu %zu) %e != %e "
			 "!= %e != %e != %e != %e\n" ,  
			 i , j , k , l ,
			 +K1xy.L1[i][j][k][l] , +K1yx.L1[i][k][j][l] ,
			 -K2xy.L1[i][l][k][j] , -K2yx.L1[i][k][l][j] ,
			 -K3xy.L1[i][l][j][k] , -K3yx.L1[i][j][l][k] ) ; 
	      }
	      // L2
	      if( fabs( K1xy.L2[i][j][k][l] - K1yx.L2[i][k][j][l] ) > 1E-12 ||
		  fabs( K1xy.L2[i][j][k][l] + K2xy.L2[i][l][k][j] ) > 1E-12 ||
		  fabs( K1xy.L2[i][j][k][l] + K2yx.L2[i][k][l][j] ) > 1E-12 ||
		  fabs( K1xy.L2[i][j][k][l] + K3xy.L2[i][l][j][k] ) > 1E-12 ||
		  fabs( K1xy.L2[i][j][k][l] + K3yx.L2[i][j][l][k] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY0 L2 failure (%zu %zu %zu %zu) %e != %e "
			 "!= %e != %e != %e != %e\n" ,  
			 i , j , k , l ,
			 +K1xy.L2[i][j][k][l] , +K1yx.L2[i][k][j][l] ,
			 -K2xy.L2[i][l][k][j] , -K2yx.L2[i][k][l][j] ,
			 -K3xy.L2[i][l][j][k] , -K3yx.L2[i][j][l][k] ) ; 
	      }
	      // L3
	      if( fabs( K1xy.L3[i][j][k][l] - K1yx.L3[i][k][j][l] ) > 1E-12 ||
		  fabs( K1xy.L3[i][j][k][l] + K2xy.L3[i][l][k][j] ) > 1E-12 ||
		  fabs( K1xy.L3[i][j][k][l] + K2yx.L3[i][k][l][j] ) > 1E-12 ||
		  fabs( K1xy.L3[i][j][k][l] + K3xy.L3[i][l][j][k] ) > 1E-12 ||
		  fabs( K1xy.L3[i][j][k][l] + K3yx.L3[i][j][l][k] ) > 1E-12 ) {
		fprintf( stderr , "SYMXY0 L3 failure (%zu %zu %zu %zu) %e != %e "
			 "!= %e != %e != %e != %e\n" ,  
			 i , j , k , l ,
			 +K1xy.L3[i][j][k][l] , +K1yx.L3[i][k][j][l] ,
			 -K2xy.L3[i][l][k][j] , -K2yx.L3[i][k][l][j] ,
			 -K3xy.L3[i][l][j][k] , -K3yx.L3[i][j][l][k] ) ; 
	      }
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  return ;
}

// test that L(x-y,0) == -L(y-x,0) && L(0,x-y) == -L(0,y-x) 
__host__
static void
example15( const struct QED_kernel_temps t )
{
  const double zero[4] = { 0 , 0 , 0 , 0 } ;

  printbar() ;
  fprintf( stdout , "Testing L(x-y,0) == -L(y-x,0) && "
	   "L(0,x-y) == L(0,y-x)\n" ) ;

  start_timer() ;
  double k1[6][4][4][4] , k2[6][4][4][4] , k3[6][4][4][4] , k4[6][4][4][4] ;
  size_t X , Y ;
  double x[4] , y[4] ;
  for( X = 0 ; X < LVOLUME ; X++ ) {
    HLBL_crds( y , X ) ;
    for( Y = 0 ; Y < LVOLUME ; Y++ ) {
      HLBL_crds( x , Y ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;
      const double ymx[4] = { -xmy[0] , -xmy[1] , -xmy[2] , -xmy[3] } ;
    
      cu_pt_QED_kernel_L0( xmy , zero , t , k1 ) ;
      cu_pt_QED_kernel_L0( ymx , zero , t , k2 ) ;
      
      cu_pt_QED_kernel_L0( zero , xmy , t , k3 ) ;
      cu_pt_QED_kernel_L0( zero , ymx , t , k4 ) ;
    
      size_t i , j , k , l ;
      for( i = 0 ; i < 6 ; i++ ) {
	for( j = 0 ; j < 4 ; j++ ) {
	  for( k = 0 ; k < 4 ; k++ ) {
	    for( l = 0 ; l < 4 ; l++ ) {
	      if( fabs( k1[i][j][k][l] + k2[i][j][k][l] ) > 1E-12 ) {
		fprintf( stdout , "TEST (%zu,%zu,%zu,%zu) %e != %e \n" ,
			 i , j , k , l , k1[i][j][k][l] , -k2[i][j][k][l] ) ;
	      }
	      if( fabs( k3[i][j][k][l] + k4[i][j][k][l] ) > 1E-12 ) {
		fprintf( stdout , "TEST (%zu,%zu,%zu,%zu) %e != %e \n" ,
			 i , j , k , l , k3[i][j][k][l] , -k4[i][j][k][l] ) ;
	      }
	    }
	  }
	}
      }
      // x and y loop
    }
  }
  printbar() ;
  print_time() ;
  
  return ;
}

// test our SYMXY0 compute all kernels code again by computing
// the individual kernels and then symmetrising them
__host__
static void
example16( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking again compute_all_kernels_SYMXY0 routine\n" ) ;

  start_timer() ;
  size_t X , Y ;
  double x[4] , y[4] ;
  for( X = 0 ; X < LVOLUME ; X++ ) {
    HLBL_crds( y , X ) ;
    for( Y = 0 ; Y < LVOLUME ; Y++ ) {
      HLBL_crds( x , Y ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;
      const double ymx[4] = { -xmy[0] , -xmy[1] ,
			      -xmy[2] , -xmy[3] } ;

      struct QED_Kernels K1 , K2 , K3 ;
      struct Kernels K ;
      cu_pt_compute_all_kernels( x , y , t , &K1 ) ;
      cu_pt_compute_all_kernels( x , xmy , t , &K2 ) ;
      cu_pt_compute_all_kernels( y , ymx , t , &K3 ) ;
      cu_pt_compute_all_kernels_SYMXY0_v2( x , y , t , &K ) ;
      
      size_t i , j , k , l ;
      double L = 0 ;
      for( i = 0 ; i < 6 ; i++ ) {
	for( j = 0 ; j < 4 ; j++ ) {
	  for( k = 0 ; k < 4 ; k++ ) {
	    for( l = 0 ; l < 4 ; l++ ) {
	      // test the L0
	      L = ( +K1.L0.xy[i][j][k][l] + K1.L0.yx[i][k][j][l]
		    -K2.L0.xy[i][l][k][j] - K2.L0.yx[i][k][l][j]
		    -K3.L0.xy[i][l][j][k] - K3.L0.yx[i][j][l][k] )/6. ;
	      if( fabs( L - K.L0[i][j][k][l] ) > 1E-12 ) {
		printf( "Failed L0 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , L , K.L0[i][j][k][l] ) ;
	      }
	      // test the L1
	      L = ( +K1.L1.xy[i][j][k][l] + K1.L1.yx[i][k][j][l]
		    -K2.L1.xy[i][l][k][j] - K2.L1.yx[i][k][l][j]
		    -K3.L1.xy[i][l][j][k] - K3.L1.yx[i][j][l][k] )/6. ;
	      if( fabs( L - K.L1[i][j][k][l] ) > 1E-12 ) {
		printf( "Failed L1 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , L , K.L1[i][j][k][l] ) ;
	      }
	      // test the L2
	      L = ( +K1.L2.xy[i][j][k][l] + K1.L2.yx[i][k][j][l]
		    -K2.L2.xy[i][l][k][j] - K2.L2.yx[i][k][l][j]
		    -K3.L2.xy[i][l][j][k] - K3.L2.yx[i][j][l][k] )/6. ;
	      if( fabs( L - K.L2[i][j][k][l] ) > 1E-12 ) {
		printf( "Failed L2 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , L , K.L2[i][j][k][l] ) ;
	      }
	      // test the L3
	      L = ( +K1.L3.xy[i][j][k][l] + K1.L3.yx[i][k][j][l]
		    -K2.L3.xy[i][l][k][j] - K2.L3.yx[i][k][l][j]
		    -K3.L3.xy[i][l][j][k] - K3.L3.yx[i][j][l][k] )/6. ;
	      if( fabs( L - K.L3[i][j][k][l] ) > 1E-12 ) {
		printf( "Failed L3 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , L , K.L3[i][j][k][l] ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  
  return ;
}

// test our SYMXY0 compute all kernels code again by computing
// the individual kernels and then symmetrising them
__host__
static void
example17( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking the con_kernels routine\n" ) ;

  start_timer() ;
  size_t X , Y ;
  double x[4] , y[4] ;
  for( X = 0 ; X < LVOLUME ; X++ ) {
    HLBL_crds( y , X ) ;
    for( Y = 0 ; Y < LVOLUME ; Y++ ) {
      HLBL_crds( x , Y ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;

      struct QED_Kernels K1 , K2 , K3 ;
      cu_pt_compute_con_kernels( x , y , t , &K1 ) ;

      cu_pt_compute_all_kernels( x , y , t , &K2 ) ;
      cu_pt_compute_all_kernels( x , xmy , t , &K3 ) ;
      
      size_t i , j , k , l ;
      double L = 0 ;
      for( i = 0 ; i < 6 ; i++ ) {
	for( j = 0 ; j < 4 ; j++ ) {
	  for( k = 0 ; k < 4 ; k++ ) {
	    for( l = 0 ; l < 4 ; l++ ) {
	      L = K2.L0.xy[i][j][k][l]+K2.L0.yx[i][k][j][l]
		-K3.L0.xy[i][l][k][j] ;
	      // test the L0
	      if( fabs( K1.L0.xy[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L0 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l ,	K1.L0.xy[i][j][k][l] , L ) ;
	      }
	      // test the L1
	      L = K2.L1.xy[i][j][k][l]+K2.L1.yx[i][k][j][l]
		-K3.L1.xy[i][l][k][j] ;
	      if( fabs( K1.L1.xy[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L1 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , K1.L1.xy[i][j][k][l] , L ) ;
	      }
	      // test the L2
	      L = K2.L2.xy[i][j][k][l]+K2.L2.yx[i][k][j][l]
		-K3.L2.xy[i][l][k][j] ;
	      if( fabs( K1.L2.xy[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L2 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , K1.L2.xy[i][j][k][l] , L ) ;
	      }
	      // test the L3
	      L = K2.L3.xy[i][j][k][l]+K2.L3.yx[i][k][j][l]
		-K3.L3.xy[i][l][k][j] ;
	      if( fabs( K1.L3.xy[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L3 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l ,	K1.L3.xy[i][j][k][l] , L ) ;
	      }
	      // test the L0
	      L = K3.L0.xy[i][l][k][j] ;
	      if( fabs( K1.L0.yx[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L0 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , K1.L0.yx[i][j][k][l] , L ) ;
	      }
	      // test the L1
	      L = K3.L1.xy[i][l][k][j] ;
	      if( fabs( K1.L1.yx[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L1 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l ,	K1.L1.yx[i][j][k][l] , L ) ;
	      }
	      // test the L2
	      L = K3.L2.xy[i][l][k][j] ;
	      if( fabs( K1.L2.yx[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L2 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l , K1.L2.yx[i][j][k][l] , L ) ;
	      }
	      // test the L3
	      L = K3.L3.xy[i][l][k][j] ;
	      if( fabs( K1.L3.yx[i][j][k][l] - L ) > 1E-12 ) {
		printf( "Failed L3 (%zu,%zu,%zu,%zu) %e %e\n" ,
			i , j , k , l ,	K1.L3.yx[i][j][k][l] , L ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  
  return ;
}

// test the con_kernels_v2 routines
__host__
static void
example18( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking the con_kernels routine\n" ) ;


  start_timer() ;
  size_t X , Y ;
  double x[4] , y[4] ;
  for( X = 0 ; X < LVOLUME ; X++ ) {
    HLBL_crds( y , X ) ;
    for( Y = 0 ; Y < LVOLUME ; Y++ ) {
      HLBL_crds( x , Y ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;

      struct Kernels K1[3] ;
      struct QED_Kernels K2 , K3 ;
      cu_pt_compute_con_kernels_v2( x , y , t , K1 ) ;

      cu_pt_compute_all_kernels( x , y , t , &K2 ) ;
      cu_pt_compute_all_kernels( x , xmy , t , &K3 ) ;
      
      size_t i , j , k , l ;
      double L = 0 ;
      for( i = 0 ; i < 6 ; i++ ) {
	for( j = 0 ; j < 4 ; j++ ) {
	  for( k = 0 ; k < 4 ; k++ ) {
	    for( l = 0 ; l < 4 ; l++ ) {
	      // L0 kernel
	      if( fabs( K2.L0.xy[i][j][k][l] - K1[0].L0[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L0 A test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K2.L0.xy[i][j][k][l] , K1[0].L0[i][j][k][l] ) ;
	      }
	      if( fabs( K3.L0.xy[i][l][k][j] - K1[1].L0[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L0 B test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K3.L0.xy[i][l][k][j] , K1[2].L0[i][j][k][l] ) ;
	      }
	      L = K2.L0.xy[i][j][k][l] + K2.L0.yx[i][k][j][l] - K3.L0.xy[i][l][k][j] ;
	      if( fabs( L - K1[2].L0[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L0 C test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , L , K1[2].L0[i][j][k][l] ) ;
	      }

	      // L1 kernel
	      if( fabs( K2.L1.xy[i][j][k][l] - K1[0].L1[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L1 A test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K2.L1.xy[i][j][k][l] , K1[0].L1[i][j][k][l] ) ;
	      }
	      if( fabs( K3.L1.xy[i][l][k][j] - K1[1].L1[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L1 B test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K3.L1.xy[i][l][k][j] , K1[1].L1[i][j][k][l] ) ;
	      }
	      L = K2.L1.xy[i][j][k][l] + K2.L1.yx[i][k][j][l] - K3.L1.xy[i][l][k][j] ;
	      if( fabs( L - K1[2].L1[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L1 C test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , L , K1[2].L1[i][j][k][l] ) ;
	      }

	      // L2 kernel
	      if( fabs( K2.L2.xy[i][j][k][l] - K1[0].L2[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L2 A test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K2.L2.xy[i][j][k][l] , K1[0].L2[i][j][k][l] ) ;
	      }
	      if( fabs( K3.L2.xy[i][l][k][j] - K1[1].L2[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L2 B test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K3.L2.xy[i][l][k][j] , K1[1].L2[i][j][k][l] ) ;
	      }
	      L = K2.L2.xy[i][j][k][l] + K2.L2.yx[i][k][j][l] - K3.L2.xy[i][l][k][j] ;
	      if( fabs( L - K1[2].L2[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L2 C test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , L , K1[2].L2[i][j][k][l] ) ;
	      }

	      // L3 kernel
	      if( fabs( K2.L2.xy[i][j][k][l] - K1[0].L2[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L3 A test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K2.L3.xy[i][j][k][l] , K1[0].L3[i][j][k][l] ) ;
	      }
	      if( fabs( K3.L2.xy[i][l][k][j] - K1[1].L2[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L3 B test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , K3.L3.xy[i][l][k][j] , K1[1].L3[i][j][k][l] ) ;
	      }
	      L = K2.L3.xy[i][j][k][l] + K2.L3.yx[i][k][j][l] - K3.L3.xy[i][l][k][j] ;
	      if( fabs( L - K1[2].L3[i][j][k][l] ) > 1E-12 ) {
		fprintf( stderr , "L3 C test failed (%zu,%zu,%zu,%zu) %e %e\n" ,
			 i , j , k , l , L , K1[2].L3[i][j][k][l] ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  
  return ;
}

// test the swap_munu_Lyx routine
__host__
static void
example19( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking the swap_munu_Lyx routine\n" ) ;

  start_timer() ;
  size_t i , j , rhosig , mu , nu , lambda ;
  double x[4] , y[4] ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    HLBL_crds( y , i ) ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      HLBL_crds( x , j ) ;
      
      struct QED_Kernels K2 , K3 ;
      cu_pt_compute_all_kernels( x , y , t , &K2 ) ;
      cu_pt_compute_all_kernels( x , y , t , &K3 ) ;
      swap_munu_Lyx( &K3 ) ;
      
      double L = 0 ;
      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      // L0 test
	      L = K2.L0.yx[rhosig][mu][nu][lambda] - K3.L0.yx[rhosig][nu][mu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "L0 swap test failed (%zu,%zu,%zu,%zu) %e %e != %e\n" ,
			 rhosig , mu , nu , lambda , L ,
			 K2.L0.yx[rhosig][mu][nu][lambda] ,
			 K3.L0.yx[rhosig][nu][mu][lambda] ) ;
	      }
	      // L1 test
	      L = K2.L1.yx[rhosig][mu][nu][lambda] - K3.L1.yx[rhosig][nu][mu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "L0 swap test failed (%zu,%zu,%zu,%zu) %e %e != %e\n" ,
			 rhosig , mu , nu , lambda , L ,
			 K2.L1.yx[rhosig][mu][nu][lambda] ,
			 K3.L1.yx[rhosig][nu][mu][lambda] ) ;
	      }
	      // L2 test
	      L = K2.L2.yx[rhosig][mu][nu][lambda] - K3.L2.yx[rhosig][nu][mu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "L0 swap test failed (%zu,%zu,%zu,%zu) %e %e != %e\n" ,
			 rhosig , mu , nu , lambda , L ,
			 K2.L2.yx[rhosig][mu][nu][lambda] ,
			 K3.L2.yx[rhosig][nu][mu][lambda] ) ;
	      }
	      // L3 test
	      L = K2.L3.yx[rhosig][mu][nu][lambda] - K3.L3.yx[rhosig][nu][mu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "L0 swap test failed (%zu,%zu,%zu,%zu) %e %e != %e\n" ,
			 rhosig , mu , nu , lambda , L ,
			 K2.L3.yx[rhosig][mu][nu][lambda] ,
			 K3.L3.yx[rhosig][nu][mu][lambda] ) ;
	      }	      
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  return ;
}

// test the  routine
__host__
static void
example20( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking the compute_con_kernelsM_L2\n" ) ;

  start_timer() ;
  double x[4] , y[4] ; 
  size_t i , j , rhosig , mu , nu , lambda ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    HLBL_crds( y , i ) ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      HLBL_crds( x , j ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;
      // todo :: a sensible test here
      double M[4] = { 0. , 0.5 , 1.0 , 1.5 } ;
      struct Kernels K[3] ;
      cu_pt_compute_con_kernelsM_L2( M , x , y , t , K ) ;

      // test against L2 for each M and xy,yx,x-y combo
      double K1[6][4][4][4] , K2[6][4][4][4] , K3[6][4][4][4] ;

      cu_pt_QED_Mkernel_L2( M[0] , x , y , t , K1 ) ;
      cu_pt_QED_Mkernel_L2( M[0] , y , x , t , K2 ) ;
      cu_pt_QED_Mkernel_L2( M[0] , x , xmy , t , K3 ) ;   
      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      if( fabs( K[0].L0[rhosig][mu][nu][lambda] - K1[rhosig][mu][nu][lambda] ) > 1E-12 ) {
		fprintf( stderr , "K1 error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[0].L0[rhosig][mu][nu][lambda] ,K1[rhosig][mu][nu][lambda] ) ;
	      }
	      if( fabs( K[1].L0[rhosig][mu][nu][lambda] - K3[rhosig][lambda][nu][mu] ) > 1E-12 ) {
		fprintf( stderr , "K3 error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[1].L0[rhosig][mu][nu][lambda] , K3[rhosig][lambda][nu][mu] ) ;
	      }
	      const double sum = K1[rhosig][mu][nu][lambda] + K2[rhosig][nu][mu][lambda] - K3[rhosig][lambda][nu][mu] ;
	      if( fabs( K[2].L0[rhosig][mu][nu][lambda] - sum ) > 1E-12 ) {
		fprintf( stderr , "Sum error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[2].L0[rhosig][mu][nu][lambda] , sum ) ;
	      }
	    }
	  }
	}
      }

      //
      cu_pt_QED_Mkernel_L2( M[1] , x , y , t , K1 ) ;
      cu_pt_QED_Mkernel_L2( M[1] , y , x , t , K2 ) ;
      cu_pt_QED_Mkernel_L2( M[1] , x , xmy , t , K3 ) ;
      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      if( fabs( K[0].L1[rhosig][mu][nu][lambda] - K1[rhosig][mu][nu][lambda] ) > 1E-12 ) {
		fprintf( stderr , "K1 error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[0].L1[rhosig][mu][nu][lambda] , K1[rhosig][mu][nu][lambda] ) ;
	      }
	      if( fabs( K[1].L1[rhosig][mu][nu][lambda] - K3[rhosig][lambda][nu][mu] ) > 1E-12 ) {
		fprintf( stderr , "K3 error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[1].L1[rhosig][mu][nu][lambda] , K3[rhosig][lambda][nu][mu] ) ;
	      }
	      const double sum = K1[rhosig][mu][nu][lambda] + K2[rhosig][nu][mu][lambda] - K3[rhosig][lambda][nu][mu] ;
	      if( fabs( K[2].L1[rhosig][mu][nu][lambda] - sum ) > 1E-12 ) {
		fprintf( stderr , "Sum error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[2].L1[rhosig][mu][nu][lambda] , sum ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  print_time() ;
  return ;
}

// test the compute_all_Mkernels
__host__
static void
example21( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking the compute_all_Mkernels routine\n" ) ;

  start_timer() ;
  double x[4] , y[4] ; 
  size_t i , j , rhosig , mu , nu , lambda ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    HLBL_crds( y , i ) ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      HLBL_crds( x , j ) ;

      const double M[4] = { 0.0 , 0.4 , 0.8 , 1.2 } ;
      struct QED_Kernels K ;
      cu_pt_compute_all_Mkernels( M , x , y , t , &K ) ;

      double L1[6][4][4][4] KQED_ALIGN , L2[6][4][4][4] KQED_ALIGN ;
      cu_pt_QED_kernel_L2( x , y , t , L1 ) ;
      cu_pt_QED_kernel_L2( y , x , t , L2 ) ;

      double L = 0.0 ;
      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      L = L1[rhosig][mu][nu][lambda] - K.L0.xy[rhosig][mu][nu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "Compute all Mkernels failed L0.xy (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , L1[rhosig][mu][nu][lambda] ,
			 K.L0.xy[rhosig][mu][nu][lambda] ) ;
	      }
	      L = L2[rhosig][mu][nu][lambda] - K.L0.yx[rhosig][mu][nu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "Compute all Mkernels failed L0.yx (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , L2[rhosig][mu][nu][lambda] ,
			 K.L0.yx[rhosig][mu][nu][lambda] ) ;
	      }
	    }
	  }
	}
      }
      
      // M != 0
      cu_pt_QED_Mkernel_L2( M[1] , x , y , t , L1 ) ;
      cu_pt_QED_Mkernel_L2( M[1] , y , x , t , L2 ) ;
      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      L = L1[rhosig][mu][nu][lambda] - K.L1.xy[rhosig][mu][nu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "Compute all Mkernels failed L0.xy (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , L1[rhosig][mu][nu][lambda] ,
			 K.L1.xy[rhosig][mu][nu][lambda] ) ;
	      }
	      L = L2[rhosig][mu][nu][lambda] - K.L1.yx[rhosig][mu][nu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "Compute all Mkernels failed L0.yx (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , L2[rhosig][mu][nu][lambda] ,
			 K.L1.yx[rhosig][mu][nu][lambda] ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  printbar() ;
  print_time() ;
  return ;
}

// test the compute_all_Mkernels
__host__
static void
example22( const struct QED_kernel_temps t )
{
  printbar() ;
  start_timer() ;
  fprintf( stdout , "Testing for thread safety\n" ) ;
#if (defined HAVE_OMP_H)
  const double x[4] = { 1 , 2 , 3 , 4 }, y[4] = { 1 , 2 , 2 , 2 } ; 
#pragma omp parallel
  {
    char str[64] ;
    sprintf( str , "Thread_%d" , omp_get_thread_num() ) ;
    FILE *file = fopen( str , "w" ) ;
    double K1[6][4][4][4] KQED_ALIGN ;
    double K2[6][4][4][4] KQED_ALIGN ;
    cu_pt_QED_kernel_L0( x , y , t , K1 ) ;
    cu_pt_QED_kernel_L0( y , x , t , K2 ) ;

    print_kernels( (const double*)K1 , (const double*)K2 , file ) ;
    
    fclose( file ) ;
  }
#endif
  print_time() ;
  return ;
}

__host__
static void
example23( const struct QED_kernel_temps t )
{  
  printbar() ;
  fprintf( stdout , "Checking the compute_all_Mkernels_v2 routine\n" ) ;

  start_timer() ;
  double x[4] , y[4] ; 
  size_t i , j , rhosig , mu , nu , lambda ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    HLBL_crds( y , i ) ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      HLBL_crds( x , j ) ;

      const double M[4] = { 0.0 , 0.4 , 0.8 , 1.2 } ;
      struct QED_Kernels K1 , K2 ;
      cu_pt_compute_all_Mkernels( M , x , y , t , &K1 ) ;
      cu_pt_compute_all_Mkernels_v2( M , x , y , t , &K2 ) ;

      double L1[6][4][4][4] KQED_ALIGN ;
      const double my[4] = { -y[0] , -y[1] , -y[2] , -y[3] } ;
      cu_pt_QED_Mkernel_L2( M[1] , x , my , t , L1 ) ;

      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {

	      // test the first term
	      double L =
		K1.L1.xy[rhosig][mu][nu][lambda] +
		K1.L1.yx[rhosig][nu][mu][lambda] -
		K2.L1.xy[rhosig][mu][nu][lambda] ;

	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "Compute all Mkernels v2 failed L1.xy (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda ,
			 K1.L1.xy[rhosig][mu][nu][lambda] +
			 K1.L1.yx[rhosig][nu][mu][lambda] ,
			 K2.L1.xy[rhosig][mu][nu][lambda] ) ;
	      }

	      L = L1[rhosig][mu][nu][lambda] - K2.L1.yx[rhosig][mu][nu][lambda] ;
	      if( fabs( L ) > 1E-12 ) {
		fprintf( stderr , "Compute all Mkernels v2 failed L1.yx (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda ,
			 L1[rhosig][mu][nu][lambda] ,
			 K2.L1.yx[rhosig][mu][nu][lambda] ) ;
	      }
	      
	    }
	  }
	}
      }
      //
    }
  }
  print_time() ;
  
  return ;
}

// test the routine
__host__
static void
example24( const struct QED_kernel_temps t )
{
  printbar() ;
  fprintf( stdout , "Checking the compute_sub_kernelsM_L2\n" ) ;

  start_timer() ;
  double x[4] , y[4] ; 
  size_t i , j , rhosig , mu , nu , lambda ;
  for( i = 0 ; i < LVOLUME ; i++ ) {
    HLBL_crds( y , i ) ;
    for( j = 0 ; j < LVOLUME ; j++ ) {
      HLBL_crds( x , j ) ;

      const double xmy[4] = { x[0]-y[0] , x[1]-y[1] ,
			      x[2]-y[2] , x[3]-y[3] } ;
      const double M[4] = { 0. , 0.5 , 1.0 , 1.5 } ;
      struct Kernels K[3] ;
      cu_pt_compute_sub_kernelsM_L2( M , x , y , t , K ) ;

      // test against L2 for each M and xy,yx,x-y combo
      double K1[6][4][4][4] , K2[6][4][4][4] , K3[6][4][4][4] ;
      cu_pt_QED_Mkernel_L2( M[0] , x , y , t , K1 ) ;
      cu_pt_QED_Mkernel_L2( M[0] , y , x , t , K2 ) ;
      cu_pt_QED_Mkernel_L2( M[0] , x , xmy , t , K3 ) ;

      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      if( fabs( K[0].L0[rhosig][mu][nu][lambda] - K1[rhosig][mu][nu][lambda] ) > 1E-12 ) {
		fprintf( stderr , "K1 error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda ,
			 K[0].L0[rhosig][mu][nu][lambda] ,K1[rhosig][mu][nu][lambda] ) ;
	      }
	      if( fabs( K[1].L0[rhosig][mu][nu][lambda] - K3[rhosig][nu][lambda][mu] ) > 1E-12 ) {
		fprintf( stderr , "K3 error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda ,
			 K[1].L0[rhosig][mu][nu][lambda] , K3[rhosig][nu][lambda][mu] ) ;
	      }
	      const double sum = K1[rhosig][mu][nu][lambda] + K2[rhosig][nu][mu][lambda] - K3[rhosig][nu][lambda][mu] ;
	      if( fabs( K[2].L0[rhosig][mu][nu][lambda] - sum ) > 1E-12 ) {
		fprintf( stderr , "Sum error M[0] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[2].L0[rhosig][mu][nu][lambda] , sum ) ;
	      }
	    }
	  }
	}
      }

      //
      cu_pt_QED_Mkernel_L2( M[1] , x , y , t , K1 ) ;
      cu_pt_QED_Mkernel_L2( M[1] , y , x , t , K2 ) ;
      cu_pt_QED_Mkernel_L2( M[1] , x , xmy , t , K3 ) ;
      for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
	for( mu = 0 ; mu < 4 ; mu++ ) {
	  for( nu = 0 ; nu < 4 ; nu++ ) {
	    for( lambda = 0 ; lambda < 4 ; lambda++ ) {
	      if( fabs( K[0].L1[rhosig][mu][nu][lambda] - K1[rhosig][mu][nu][lambda] ) > 1E-12 ) {
		fprintf( stderr , "K1 error M[1] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda ,
			 K[0].L1[rhosig][mu][nu][lambda] , K1[rhosig][mu][nu][lambda] ) ;
	      }
	      if( fabs( K[1].L1[rhosig][mu][nu][lambda] - K3[rhosig][nu][lambda][mu] ) > 1E-12 ) {
		fprintf( stderr , "K3 error M[1] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda ,
			 K[1].L1[rhosig][mu][nu][lambda] , K3[rhosig][nu][lambda][mu] ) ;
	      }
	      const double sum = K1[rhosig][mu][nu][lambda] + K2[rhosig][nu][mu][lambda] - K3[rhosig][nu][lambda][mu] ;
	      if( fabs( K[2].L1[rhosig][mu][nu][lambda] - sum ) > 1E-12 ) {
		fprintf( stderr , "Sum error M[1] (%zu,%zu,%zu,%zu) %e != %e\n" ,
			 rhosig , mu , nu , lambda , K[2].L0[rhosig][mu][nu][lambda] , sum ) ;
	      }
	    }
	  }
	}
      }
    }
  }
  print_time() ;
  return ;
}

// gives us an idea of how long it takes for V full kernel calls
__host__
static void
stress_test( const struct QED_kernel_temps t )
{  
  const size_t V = 8*8*8*16 ;
  size_t i ;
  const double x[4] = { 1 , 2, 3, 4 } ;
  const double y[4] = { 1 , 1, 1, 1 };

  printbar() ;
  fprintf( stdout , "Timing test for %zu calls to all kernels\n" , V ) ;
  start_timer() ;
  
  const double M[4] = { 0. , 0.5 , 1.0 , 1.5 } ;
  for( i = 0 ; i < V ; i++ ) {
    struct Kernels K[3] ;
    cu_pt_compute_con_kernelsM_L2( M , x , y , t , K ) ;
  }
  print_time() ;
#if (defined HAVE_OMP_H)
  printbar() ;
  fprintf( stdout , "Multithreaded timing with %d threads\n" , omp_get_max_threads() ) ;
  start_timer() ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < V ; i++ ) {
    struct Kernels K[3] ;
    cu_pt_compute_con_kernelsM_L2( M , x , y , t , K ) ;
  }
  print_time() ;
#endif

  printbar() ;
  return ;
}

// just an example of how to call the code and a few simple examples/tests
__host__
int
main( void )
{  
  // initialise G factors
  struct QED_kernel_temps t ;

  start_timer() ;

  if( initialise( &t ) ) {
    goto memfree ;
  }
  
  print_time() ;

  print_header( 1 ) ;
  example1( t ) ;

  print_header( 2 ) ;
  example2( t ) ;

  print_header( 3 ) ;
  example3( t ) ;

  print_header( 4 ) ;
  example4( t ) ;

  print_header( 5 ) ;
  example5( t ) ;

  print_header( 6 ) ;
  example6( t ) ;

  print_header( 7 ) ;
  example7( t ) ;
  
  print_header( 9 ) ;
  example9( t ) ;

  print_header( 10 ) ;
  example10( t ) ;

  print_header( 11 ) ;
  example11( t ) ;

  print_header( 12 ) ;
  example12( t ) ;

  print_header( 13 ) ;
  example13( t ) ;

  print_header( 14 ) ;
  example14( t ) ;

  print_header( 15 ) ;
  example15( t ) ;

  print_header( 16 ) ;
  example16( t ) ;

  print_header( 17 ) ;
  example17( t ) ;

  print_header( 18 ) ;
  example18( t ) ;

  print_header( 19 ) ;
  example19( t ) ;

  print_header( 20 ) ;
  example20( t ) ;

  print_header( 21 ) ;
  example21( t ) ;

  // SKIP: Testing thread safety
  // print_header( 22 ) ;
  // example22( t ) ;

  print_header( 23 ) ;
  example23( t ) ;
  
  print_header( 24 ) ;
  example24( t ) ;

  // SKIP: Need to reformulate appropriately for GPU to stress test
  // printbar() ;
  // fprintf( stdout , "=== STRESS TEST ===\n" ) ;
  // stress_test( t ) ;

 memfree :
  
  free_QED_temps( &t ) ;

  fprintf( stdout , "KQED| Freed QED kernel temporaries\n" ) ;
  
  return 0 ;
}
