/**
   @file io.c
   @brief reading and writing stuff
 */
#include "KQED.h"

#include "crc32c.h"      // DML_checksum_accum
#include "GLU_bswap.h"   // byte swaps if we are different endian to file format

#define FILE_FF "/PRECOMP/FFxy_single_cksum.bin"
#define FILE_TAYLORX "/PRECOMP/taylorx_cksum.bin"
#define FILE_TAYLORY "/PRECOMP/taylory_cksum.bin" 

__host__
static int
FREAD32( void *p , const size_t size , const size_t length , FILE *f )
{
  if( fread( p , size , length , f ) != length ) {
    fprintf( stderr , "FREAD32 failure\n" ) ;
    return 1 ;
  }
#ifdef WORDS_BIGENDIAN
  bswap_32( length , p ) ;
#endif
  return 0 ;
}

__host__
static int
FREAD64( void *p , const size_t size , const size_t length , FILE *f )
{
  if( fread( p , size , length , f ) != length ) {
    fprintf( stderr , "FREAD32 failure\n" ) ;
    return 1 ;
  }
#ifdef WORDS_BIGENDIAN
  bswap_32( length , p ) ;
#endif
  return 0 ;
}

// read the (single-precision) binary file allocating Grid objects
__host__
int
read_ff( struct Grid_coeffs *Grid )
{
  char filestr[ strlen( LIB_PATH ) + strlen( FILE_FF ) + 1 ] ;
  sprintf( filestr , "%s%s" , LIB_PATH , FILE_FF ) ;
  FILE *fr = fopen( filestr , "rb" ) ;

  if( fr == NULL ) {
    fprintf( stderr , "Cannot find %s\n" , filestr ) ;
    return 1 ;
  }

  int moniker ;
  FREAD32( &moniker , sizeof( int ) , 1 , fr ) ;

  if( moniker != 816968 ) {
    fprintf( stderr , "|IO| missing magic number in FFxy_single.bin\n" ) ;
    return 1 ;
  }

  // XX
  FREAD32( &Grid -> nstpx , sizeof( int ) , 1 , fr ) ;
  double *XX = (double*)malloc( Grid -> nstpx * sizeof( double ) ) ;
  FREAD64( XX , sizeof( double ) , Grid -> nstpx , fr ) ;
  uint32_t cksumXX[2] = { 0 , 0 } , cksumXX_r[2] ;
  DML_checksum_accum_crc32c( &cksumXX[0] , &cksumXX[1] , 
			     0 , (char*)Grid -> XX ,
			     Grid -> nstpx*sizeof(double) ) ;
  FREAD32( cksumXX_r , sizeof( uint32_t ) , 2 , fr ) ;
  if( cksumXX[0] != cksumXX_r[0] || cksumXX[1] != cksumXX_r[1] ) {
    fprintf( stderr , "Computed and read XX checksums do not match!" ) ; 
    fprintf( stderr , "XX_a %u != %u\n" , cksumXX[0] , cksumXX_r[0] ) ;
    fprintf( stderr , "XX_b %u != %u\n" , cksumXX[1] , cksumXX_r[1] ) ;
    return 1 ;
  }
  checkCudaErrors(cudaMalloc(&Grid->XX, Grid -> nstpx * sizeof( double )));
  checkCudaErrors(cudaMemcpy(
      Grid->XX, XX, Grid -> nstpx * sizeof( double ), cudaMemcpyHostToDevice));
  free(XX);

  // allocate nfx
  checkCudaErrors(cudaMalloc(&Grid->nfx, Grid -> nstpx * sizeof( int ) )) ;
  int* nfx = (int*)malloc( Grid -> nstpx * sizeof( int ) ) ;

  // YY
  FREAD32( &Grid -> nstpy , sizeof( int ) , 1 , fr ) ;
  double *YY = (double*)malloc( Grid -> nstpy * sizeof( double ) ) ;
  FREAD64( YY , sizeof( double ) , Grid -> nstpy , fr ) ;
  uint32_t cksumYY[2] = { 0 , 0 } , cksumYY_r[2] ;
  DML_checksum_accum_crc32c( &cksumYY[0] , &cksumYY[1] , 
			     0 , (char*)Grid -> YY ,
			     Grid -> nstpy*sizeof(double) ) ;
  FREAD32( cksumYY_r , sizeof( uint32_t ) , 2 , fr ) ;
  if( cksumYY[0] != cksumYY_r[0] || cksumYY[1] != cksumYY_r[1] ) {
    fprintf( stderr , "Computed and read YY checksums do not match!" ) ; 
    fprintf( stderr , "XX_a %u != %u\n" , cksumYY[0] , cksumYY_r[0] ) ;
    fprintf( stderr , "XX_b %u != %u\n" , cksumYY[1] , cksumYY_r[1] ) ;
    return 1 ;
  }
  checkCudaErrors(cudaMalloc(&Grid->YY, Grid -> nstpy * sizeof( double )));
  checkCudaErrors(cudaMemcpy(
      Grid->YY, YY, Grid -> nstpy * sizeof( double ), cudaMemcpyHostToDevice));
  free(YY);

#ifdef VERBOSE
  fprintf( stdout , "Nstpx %d :: Nstpy %d\n" , Grid -> nstpx , Grid -> nstpy ) ;
#endif
  
  Grid -> xstp = ( Grid -> XX[1] - Grid -> XX[0] ) ;
  Grid -> ystp = ( Grid -> YY[1] - Grid -> YY[0] ) ;    

#ifdef VERBOSE
  fprintf( stdout , "xstp %e :: ystp %e\n" , Grid -> xstp , Grid -> ystp ) ;
#endif
  
  FREAD32( &Grid -> Nffa , sizeof( int ) , 1 , fr ) ;
  Grid -> nfx_max = 0;

  float**** Ffm = (float****)malloc( Grid -> Nffa*sizeof(float***)) ;
  float**** Ffp = (float****)malloc( Grid -> Nffa*sizeof(float***)) ;
  
  int i , j , k , NX , NY , nx ;
  size_t rank = 0 ;
  uint32_t cksumFfm_a = 0 , cksumFfm_b = 0 ;
  uint32_t cksumFfp_a = 0 , cksumFfp_b = 0 ;
  for( i = 0 ; i < Grid -> Nffa ; i++ ) {
    FREAD32( &NX , sizeof( int ) , 1 , fr ) ;

    Ffm[i] = (float***)malloc( Grid -> nstpx*sizeof(float**)) ;
    Ffp[i] = (float***)malloc( Grid -> nstpx*sizeof(float**)) ; 

    if( NX != Grid -> nstpx ) {
      fprintf( stderr , "|IO| file is funny nstpx\n" ) ;
      return 1 ;
    }
    
    for( j = 0 ; j < NX ; j++ ) {
      FREAD32( &NY , sizeof( uint32_t ) , 1 , fr ) ;
      
      if( NY != Grid -> nstpy ) {
	fprintf( stderr , "|IO| file is funny nstpy\n" ) ;
	return 1 ;
      }

      Ffm[i][j] = (float**)malloc( Grid-> nstpy*sizeof(float*)) ;
      Ffp[i][j] = (float**)malloc( Grid-> nstpy*sizeof(float*)) ;
      
      for( k = 0 ; k < NY ; k++ ) {
	FREAD32( &nx , sizeof( int ) , 1 , fr ) ;

	Ffm[i][j][k] = (float*)malloc(nx*sizeof(float)) ;
	Ffp[i][j][k] = (float*)malloc(nx*sizeof(float)) ;

	// set the length of nfx
        if (k > 0 && nx != nfx[j]) {
          fprintf( stderr, "nx mismatch\n" );
          return 1 ;
        }
	nfx[j] = nx ;
        if (nx > Grid -> nfx_max) {
          Grid -> nfx_max = nx;
        }
	
	FREAD32( Ffm[i][j][k] , sizeof( float ) , nx , fr ) ;
	DML_checksum_accum_crc32c( &cksumFfm_a , &cksumFfm_b , 
				   rank , (char*) Ffm[i][j][k] ,
				   nx*sizeof(float) ) ;
	
	//
	FREAD32( Ffp[i][j][k] , sizeof( float ) , nx , fr ) ;
	DML_checksum_accum_crc32c( &cksumFfp_a , &cksumFfp_b , 
				   rank , (char*) Ffp[i][j][k] ,
				   nx*sizeof(float) ) ;
	//
	rank ++ ;
      }
    }
  }

  uint32_t cksum[4] ;
  FREAD32( cksum , sizeof( uint32_t ) , 4 , fr ) ;

  fclose( fr ) ;
  
  if( cksumFfm_a != cksum[0] || cksumFfm_b != cksum[1] ||
      cksumFfp_a != cksum[2] || cksumFfp_b != cksum[3] ) {
    fprintf( stderr , "Computed and read ff checksums do not match!" ) ; 
    fprintf( stderr , "Ffm_a %u != %u\n" , cksum[0] , cksumFfm_a ) ;
    fprintf( stderr , "Ffm_b %u != %u\n" , cksum[1] , cksumFfm_a ) ;
    fprintf( stderr , "Ffp_a %u != %u\n" , cksum[2] , cksumFfm_a ) ;
    fprintf( stderr , "Ffp_b %u != %u\n" , cksum[3] , cksumFfm_a ) ;
    return 1 ;
  }


  // build flattened rectangular array
  size_t sizeof_Ff = (
      Grid -> Nffa * Grid -> nstpx * Grid -> nstpy *
      Grid -> nfx_max * sizeof(float) );
  float* Ffp_rect = (float*)malloc( sizeof_Ff );
  float* Ffm_rect = (float*)malloc( sizeof_Ff );
  for( i = 0; i < Grid -> Nffa ; i++ ) {
    for( j = 0; j < Grid -> nstpx ; j++ ) {
      for( k = 0; k < Grid -> nstpy ; k++ ) {
        size_t ind = ((i * Grid -> nstpx + j) * Grid -> nstpy + k) * Grid -> nfx_max ;
        memcpy( &Ffp_rect[ind], Ffp[i][j][k], Grid -> nfx_max * sizeof(float) ) ;
        memcpy( &Ffm_rect[ind], Ffm[i][j][k], Grid -> nfx_max * sizeof(float) ) ;
        free( Ffp[i][j][k] ) ;
        free( Ffm[i][j][k] ) ;
      }
      free( Ffp[i][j] ) ;
      free( Ffm[i][j] ) ;
    }
    free( Ffp[i] ) ;
    free( Ffm[i] ) ;
  }
  free( Ffp );
  free( Ffm );

  checkCudaErrors(cudaMalloc( &Grid -> Ffp, sizeof_Ff )) ;
  checkCudaErrors(cudaMalloc( &Grid -> Ffm, sizeof_Ff )) ;
  checkCudaErrors(cudaMemcpy(
      Grid -> Ffp, Ffp_rect, sizeof_Ff, cudaMemcpyHostToDevice )) ;
  checkCudaErrors(cudaMemcpy(
      Grid -> Ffm, Ffm_rect, sizeof_Ff, cudaMemcpyHostToDevice )) ;
  free( Ffp_rect );
  free( Ffm_rect );

  checkCudaErrors(cudaMemcpy(
      Grid -> nfx, nfx, Grid -> Nffa * sizeof(int), cudaMemcpyHostToDevice )) ;
  free( nfx );

  return 0 ;
}

int
read_TAYLORX( struct Grid_coeffs *Grid )
{
  char filestr[ strlen( LIB_PATH ) + strlen( FILE_TAYLORX ) + 1 ] ;
  sprintf( filestr , "%s%s" , LIB_PATH , FILE_TAYLORX ) ;

  FILE *fr = fopen( filestr , "rb" ) ;
  if( fr == NULL ) {
    fprintf( stderr , "|IO| cannot open %s\n" , filestr ) ;
    return 1 ;
  }
  int moniker ;
  FREAD32( &moniker , sizeof( int ) , 1 , fr ) ;
  if( moniker != 816968 ) {
    fprintf( stderr , "|IO| misread magic number\n" ) ;
    return 1 ;
  }
  FREAD32( &Grid -> NtayY , sizeof( int ) , 1 , fr ) ;
  if( Grid -> NtayY != TX_LEN ) {
    fprintf( stderr , "|IO| taylorx.bin misread %d vs. %d\n" ,
	     Grid -> NtayY , TX_LEN ) ;
    return 1 ;
  }
  double *TX = (double*)malloc( TX_LEN * NYTAY * sizeof( double ) ) ;
  uint32_t cksuma = 0 , cksumb = 0 ;
  size_t i ;
  for( i = 0 ; i < (size_t)Grid -> NtayY ; i++ ) {
    FREAD32( &Grid -> NY_tay , sizeof( int ) , 1 , fr ) ;
    if( Grid -> NY_tay != NYTAY ) {
      fprintf( stderr , "|IO| taylorx.bin weird NY_TAY %d vs. %d\n" ,
	       Grid -> NY_tay, NYTAY ) ;
      return 1 ;
    }
    // TX[i] = malloc( Grid -> NY_tay * sizeof( double ) ) ;
    FREAD64( &TX[i * NYTAY] , sizeof( double ) , Grid -> NY_tay , fr ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , &TX[i * NYTAY] ,
			       Grid->NY_tay * sizeof(double) ) ;
  }
  checkCudaErrors(cudaMalloc( &Grid -> TX, TX_LEN * NYTAY * sizeof( double ) )) ;
  checkCudaErrors(cudaMemcpy(
      Grid -> TX, TX, TX_LEN * NYTAY * sizeof( double ), cudaMemcpyHostToDevice )) ;
  free( TX );

  uint32_t cksum[2] ;
  FREAD32( cksum , sizeof( uint32_t ) , 2 , fr ) ;

  fclose( fr ) ;
  
  if( cksum[0] != cksuma || cksum[1] != cksumb ) {
    fprintf( stderr , "Read & Computed taylorx checksums do not match!" ) ; 
    fprintf( stderr , "Ffm_a %u != %u\n" , cksum[0] , cksuma ) ;
    fprintf( stderr , "Ffm_b %u != %u\n" , cksum[1] , cksumb ) ;
    return 1 ;
  }

  return 0 ;
}

// read the file taylory
int
read_TAYLORY( struct Grid_coeffs *Grid )
{
  char filestr[ strlen( LIB_PATH ) + strlen( FILE_TAYLORY ) + 1 ] ;
  sprintf( filestr , "%s%s" , LIB_PATH , FILE_TAYLORY ) ;

  FILE *fr = fopen( filestr , "rb" ) ;
  if( fr == NULL ) {
    fprintf( stderr , "|IO| cannot open %s\n" , filestr ) ;
    return 1 ;
  }

  int moniker ;
  FREAD32( &moniker , sizeof( int ) , 1 , fr ) ;
  if( moniker != 816968 ) {
    fprintf( stderr , "|IO| misread magic number\n" ) ;
    return 1 ;
  }
  FREAD32( &Grid -> NtayX , sizeof( int ) , 1 , fr ) ;
  if( Grid -> NtayX != 14 ) {
    fprintf( stderr , "|IO| taylory.bin misread %d vs. %d\n" ,
	     Grid -> NtayX , 14 ) ;
    return 1 ;
  }
  double *TY = (double*)malloc( TY_LEN * NXTAY * sizeof( double ) ) ;
  size_t i ;
  uint32_t cksuma = 0 , cksumb = 0 ;
  for( i = 0 ; i < (size_t)Grid -> NtayX ; i++ ) {
    FREAD32( &Grid -> NX_tay , sizeof( int ) , 1 , fr ) ;
    if( Grid -> NX_tay != NXTAY ) {
      fprintf( stderr , "|IO| taylory.bin weird NX_TAY %d vs. %d\n" ,
	       Grid -> NX_tay, NXTAY ) ;
      return 1 ;
    }
    // Grid -> TY[i] = malloc( Grid -> NX_tay * sizeof( double ) ) ;
    FREAD64( &Grid -> TY[i * NXTAY] , sizeof( double ) , Grid -> NX_tay , fr ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , &Grid -> TY[i * NXTAY] ,
			       Grid->NX_tay * sizeof(double) ) ;
  }
  checkCudaErrors(cudaMalloc( &Grid -> TY, TY_LEN * NXTAY * sizeof( double ) )) ;
  checkCudaErrors(cudaMemcpy(
      Grid -> TY, TY, TY_LEN * NXTAY * sizeof( double ), cudaMemcpyHostToDevice )) ;
  free( TY );

  uint32_t cksum[2] ;
  FREAD32( cksum , sizeof( uint32_t ) , 2 , fr ) ;

  fclose( fr ) ;
  
  if( cksum[0] != cksuma || cksum[1] != cksumb ) {
    fprintf( stderr , "Read & Computed taylory checksums do not match!" ) ; 
    fprintf( stderr , "Ffm_a %u != %u\n" , cksum[0] , cksuma ) ;
    fprintf( stderr , "Ffm_b %u != %u\n" , cksum[1] , cksumb ) ;
    return 1 ;
  }
  
  return 0 ;
}

// write the form factor
void
write_ff( const struct Grid_coeffs *Grid )
{
  // write it out in one big file
  const int moniker = 816968 ;
  FILE *fw = fopen( "./PRECOMP/FFxy_single_cksum2.bin" , "wb" ) ;
  fwrite( &moniker , sizeof( int ) , 1 , fw ) ;
  
  fwrite( &Grid -> nstpx , sizeof( int ) , 1 , fw ) ;
  double *XX = (double*)malloc( Grid -> nstpx * sizeof( double ) ) ;
  checkCudaErrors(cudaMemcpy(
      XX, Grid -> XX, Grid -> nstpx * sizeof( double ), cudaMemcpyDeviceToHost )) ;
  fwrite( XX , sizeof( double ) , Grid -> nstpx , fw ) ;

  uint32_t cksumXX[2] = { 0 , 0 } ;
  DML_checksum_accum_crc32c( &cksumXX[0] , &cksumXX[1] , 
			     0 , (char*)XX ,
			     Grid -> nstpx*sizeof(double) ) ;
  fwrite( cksumXX , sizeof( uint32_t ) , 2 , fw ) ;
  free( XX ) ;
  
  fwrite( &Grid -> nstpy , sizeof( int ) , 1 , fw ) ;
  double *YY = (double*)malloc( Grid -> nstpy * sizeof( double ) ) ;
  checkCudaErrors(cudaMemcpy(
      YY, Grid -> YY, Grid -> nstpy * sizeof( double ), cudaMemcpyDeviceToHost )) ;
  fwrite( YY , sizeof( double ) , Grid -> nstpy , fw ) ;

  uint32_t cksumYY[2] = { 0 , 0 } ;
  DML_checksum_accum_crc32c( &cksumYY[0] , &cksumYY[1] ,
			     0 , (char*)YY ,
			     Grid -> nstpy*sizeof(double) ) ;
  fwrite( cksumYY , sizeof( uint32_t ) , 2 , fw ) ;
  free( YY ) ;
  
  int Nffa = Grid -> Nffa ;
  int i , j , k ;

  uint32_t cksumFfm_a = 0 , cksumFfm_b = 0 ;
  uint32_t cksumFfp_a = 0 , cksumFfp_b = 0 ;
  size_t rank = 0 ;

  size_t sizeof_Ff = (
      Grid -> Nffa * Grid -> nstpx * Grid -> nstpy *
      Grid -> nfx_max * sizeof(float) );
  float *Ffp_rect = (float*)malloc( sizeof_Ff );
  float *Ffm_rect = (float*)malloc( sizeof_Ff );
  int *nfx = (int*)malloc( Grid -> nstpx * sizeof( int ) ) ;
  checkCudaErrors(cudaMemcpy(
      Ffp_rect, Grid -> Ffp, sizeof_Ff, cudaMemcpyDeviceToHost )) ;
  checkCudaErrors(cudaMemcpy(
      Ffm_rect, Grid -> Ffm, sizeof_Ff, cudaMemcpyDeviceToHost )) ;
  checkCudaErrors(cudaMemcpy(
      nfx, Grid -> nfx, Grid -> nstpx * sizeof( int ), cudaMemcpyDeviceToHost )) ;
  fwrite( &Nffa , sizeof( int ) , 1 , fw ) ;
  for( i = 0 ; i < Nffa ; i++ ) {
    fwrite( &Grid -> nstpx , sizeof( int ) , 1 , fw ) ;
    for( j = 0 ; j < Grid -> nstpx ; j++ ) {
      int nx = nfx[j] ;
      fwrite( &Grid -> nstpy , sizeof( int ) , 1 , fw ) ;
      for( k = 0 ; k < Grid -> nstpy ; k++ ) {
	fwrite( &nx , sizeof( int ) , 1 , fw ) ;
        float* buf;
        buf = &Ffm_rect[((i*Grid->nstpx + j)*Grid->nstpy + k) * Grid -> nfx_max];
	fwrite( buf , sizeof( float ) , nx , fw ) ;
	DML_checksum_accum_crc32c( &cksumFfm_a , &cksumFfm_b , 
				   rank , (char*)buf , nx*sizeof(float) ) ;
        buf = &Ffp_rect[((i*Grid->nstpx + j)*Grid->nstpy + k) * Grid -> nfx_max];
	fwrite( buf , sizeof( float ) , nx , fw ) ;
	DML_checksum_accum_crc32c( &cksumFfp_a , &cksumFfp_b , 
				   rank , (char*)buf , nx*sizeof(float) ) ;
	rank++ ;
      }
    }
  }
  fwrite( &cksumFfm_a , sizeof( uint32_t ) , 1 , fw ) ;
  fwrite( &cksumFfm_b , sizeof( uint32_t ) , 1 , fw ) ;
  fwrite( &cksumFfp_a , sizeof( uint32_t ) , 1 , fw ) ;
  fwrite( &cksumFfp_b , sizeof( uint32_t ) , 1 , fw ) ;

  free( Ffp_rect );
  free( Ffm_rect );
  free( nfx );

  fclose( fw ) ;
  return ;
}

void
write_TAYLORX( const struct Grid_coeffs *Grid )
{
  FILE *fw = fopen( "./PRECOMP/taylorx_cksum.bin", "wb" ) ;
  if( fw == NULL ) {
    fprintf( stderr , "|IO| cannot write to file TAYLORX\n" ) ;
    return ;
  }
  const int moniker = 816968 ;
  fwrite( &moniker , sizeof( int ) , 1 , fw ) ;
  fwrite( &Grid -> NtayY , sizeof( int ) , 1 , fw ) ;
  uint32_t cksuma = 0 , cksumb = 0 ;
  size_t i ;
  double *TX  = (double*)malloc(Grid -> NtayY * Grid -> NY_tay * sizeof( double ) ) ;
  checkCudaErrors(cudaMemcpy(
      TX, Grid -> TX, Grid -> NtayY * Grid -> NY_tay * sizeof( double ), cudaMemcpyDeviceToHost )) ;
  for( i = 0 ; i < (size_t)Grid -> NtayY ; i++ ) {
    fwrite( &Grid -> NY_tay , sizeof( int ) , 1 , fw ) ;
    fwrite( &TX[i * Grid -> NY_tay] , sizeof( double ) , Grid -> NY_tay , fw  ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , &Grid -> TX[i * Grid -> NY_tay] ,
			       Grid->NY_tay * sizeof(double) ) ;
  }
  free( TX ) ;
  fwrite( &cksuma , sizeof( uint32_t ) , 1 , fw ) ;
  fwrite( &cksumb , sizeof( uint32_t ) , 1 , fw ) ;
  fclose( fw ) ;
  return ;
}

void
write_TAYLORY( const struct Grid_coeffs *Grid )
{  
  FILE *fw = fopen( "./PRECOMP/taylory_cksum.bin", "wb" ) ;
  if( fw == NULL ) {
    fprintf( stderr , "|IO| cannot write to file TAYLORY\n" ) ;
    return ;
  }
  const int moniker = 816968 ;
  fwrite( &moniker , sizeof( int ) , 1 , fw ) ;
  fwrite( &Grid -> NtayX , sizeof( int ) , 1 , fw ) ;
  
  uint32_t cksuma = 0 , cksumb = 0 ;
  size_t i ;
  double *TY  = (double*)malloc(Grid -> NtayX * Grid -> NX_tay * sizeof( double ) ) ;
  checkCudaErrors(cudaMemcpy(
      TY, Grid -> TY, Grid -> NtayX * Grid -> NX_tay * sizeof( double ), cudaMemcpyDeviceToHost )) ;
  for( i = 0 ; i < (size_t)Grid -> NtayX ; i++ ) {
    fwrite( &Grid -> NX_tay , sizeof( int ) , 1 , fw ) ;
    fwrite( &TY[i * Grid -> NX_tay] , sizeof( double ) , Grid -> NX_tay , fw  ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , &TY[i * Grid -> NX_tay] ,
			       Grid->NX_tay * sizeof(double) ) ;
  }
  free( TY ) ;
  
  fwrite( &cksuma , sizeof( uint32_t ) , 1 , fw ) ;
  fwrite( &cksumb , sizeof( uint32_t ) , 1 , fw ) ;
  
  fclose( fw ) ;

  return ;
}

#ifdef FILE_FF
  #undef FILE_FF
#endif

#ifdef FILE_TAYLORX
  #undef FILE_TAYLORX
#endif

#ifdef FILE_TAYLORY
  #undef FILE_TAYLORY
#endif
