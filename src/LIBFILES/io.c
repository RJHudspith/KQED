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
  
  FREAD32( &Grid -> nstpx , sizeof( int ) , 1 , fr ) ;
  Grid -> XX = malloc( (size_t)Grid -> nstpx * sizeof( double ) ) ;
  FREAD64( Grid -> XX , sizeof( double ) , Grid -> nstpx , fr ) ;
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
  // allocate nfx
  Grid -> nfx = malloc( (size_t)Grid -> nstpx * sizeof( int ) ) ;
  
  FREAD32( &Grid -> nstpy , sizeof( int ) , 1 , fr ) ;
  Grid -> YY = malloc( Grid -> nstpy * sizeof( double ) ) ;
  FREAD64( Grid -> YY , sizeof( double ) , Grid -> nstpy , fr ) ;
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

#ifdef VERBOSE
  fprintf( stdout , "Nstpx %d :: Nstpy %d\n" , Grid -> nstpx , Grid -> nstpy ) ;
#endif
  
  Grid -> xstp = ( Grid -> XX[1] - Grid -> XX[0] ) ;
  Grid -> ystp = ( Grid -> YY[1] - Grid -> YY[0] ) ;    

#ifdef VERBOSE
  fprintf( stdout , "xstp %e :: ystp %e\n" , Grid -> xstp , Grid -> ystp ) ;
#endif
  
  FREAD32( &Grid -> Nffa , sizeof( int ) , 1 , fr ) ;
  
  Grid -> Ffm = (float****)malloc( Grid -> Nffa*sizeof(float***)) ;
  Grid -> Ffp = (float****)malloc( Grid -> Nffa*sizeof(float***)) ; 
  
  int i , j , k , NX , NY , nx ;
  size_t rank = 0 ;
  uint32_t cksumFfm_a = 0 , cksumFfm_b = 0 ;
  uint32_t cksumFfp_a = 0 , cksumFfp_b = 0 ;
  for( i = 0 ; i < Grid -> Nffa ; i++ ) {
    FREAD32( &NX , sizeof( int ) , 1 , fr ) ;

    Grid -> Ffm[i] = (float***)malloc( Grid -> nstpx*sizeof(float**)) ;
    Grid -> Ffp[i] = (float***)malloc( Grid -> nstpx*sizeof(float**)) ; 

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

      Grid -> Ffm[i][j] = (float**)malloc( Grid-> nstpy*sizeof(float*)) ;
      Grid -> Ffp[i][j] = (float**)malloc( Grid-> nstpy*sizeof(float*)) ;
      
      for( k = 0 ; k < NY ; k++ ) {
	FREAD32( &nx , sizeof( int ) , 1 , fr ) ;

	Grid -> Ffm[i][j][k] = (float*)malloc(nx*sizeof(float)) ;
	Grid -> Ffp[i][j][k] = (float*)malloc(nx*sizeof(float)) ;

	// set the length of nfx
	Grid -> nfx[j] = nx ;
	
	FREAD32( Grid -> Ffm[i][j][k] , sizeof( float ) , nx , fr ) ;
	DML_checksum_accum_crc32c( &cksumFfm_a , &cksumFfm_b , 
				   rank , (char*)Grid -> Ffm[i][j][k] ,
				   nx*sizeof(float) ) ;
	
	//
	FREAD32( Grid -> Ffp[i][j][k] , sizeof( float ) , nx , fr ) ;
	DML_checksum_accum_crc32c( &cksumFfp_a , &cksumFfp_b , 
				   rank , (char*)Grid -> Ffp[i][j][k] ,
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
  if( Grid -> NtayY != 23 ) {
    fprintf( stderr , "|IO| taylorx.bin misread %d vs. %d\n" ,
	     Grid -> NtayY , 23 ) ;
    return 1 ;
  }
  Grid -> TX = malloc( Grid -> NtayY * sizeof( double* ) ) ;
  uint32_t cksuma = 0 , cksumb = 0 ;
  size_t i ;
  for( i = 0 ; i < (size_t)Grid -> NtayY ; i++ ) {
    FREAD32( &Grid -> NY_tay , sizeof( int ) , 1 , fr ) ;
    if( Grid -> NY_tay != 810 ) {
      fprintf( stderr , "|IO| taylorx.bin weird NY_TAY %d vs. 810\n" ,
	       Grid -> NY_tay ) ;
      return 1 ;
    }
    Grid -> TX[i] = malloc( Grid -> NY_tay * sizeof( double ) ) ;
    FREAD64( Grid -> TX[i] , sizeof( double ) , Grid -> NY_tay , fr ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , Grid -> TX[i] ,
			       Grid->NY_tay * sizeof(double) ) ;
  }

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
  Grid -> TY = malloc( Grid -> NtayX * sizeof( double* ) ) ;
  size_t i ;
  uint32_t cksuma = 0 , cksumb = 0 ;
  for( i = 0 ; i < (size_t)Grid -> NtayX ; i++ ) {
    FREAD32( &Grid -> NX_tay , sizeof( int ) , 1 , fr ) ;
    if( Grid -> NX_tay != 100 ) {
      fprintf( stderr , "|IO| taylory.bin weird NX_TAY %d vs. 100\n" ,
	       Grid -> NX_tay ) ;
      return 1 ;
    }
    Grid -> TY[i] = malloc( Grid -> NX_tay * sizeof( double ) ) ;
    FREAD64( Grid -> TY[i] , sizeof( double ) , Grid -> NX_tay , fr ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , Grid -> TY[i] ,
			       Grid->NX_tay * sizeof(double) ) ;
  }

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
  fwrite( Grid -> XX , sizeof( double ) , Grid -> nstpx , fw ) ;

  uint32_t cksumXX[2] = { 0 , 0 } ;
  DML_checksum_accum_crc32c( &cksumXX[0] , &cksumXX[1] , 
			     0 , (char*)Grid -> XX ,
			     Grid -> nstpx*sizeof(double) ) ;
  fwrite( cksumXX , sizeof( uint32_t ) , 2 , fw ) ;
  
  fwrite( &Grid -> nstpy , sizeof( int ) , 1 , fw ) ;
  fwrite( Grid -> YY , sizeof( double ) , Grid -> nstpy , fw ) ;

  uint32_t cksumYY[2] = { 0 , 0 } ;
  DML_checksum_accum_crc32c( &cksumYY[0] , &cksumYY[1] ,
			     0 , (char*)Grid -> YY ,
			     Grid -> nstpy*sizeof(double) ) ;
  fwrite( cksumYY , sizeof( uint32_t ) , 2 , fw ) ;
  
  int Nffa = Grid -> Nffa ;
  int i , j , k ;

  uint32_t cksumFfm_a = 0 , cksumFfm_b = 0 ;
  uint32_t cksumFfp_a = 0 , cksumFfp_b = 0 ;
  size_t rank = 0 ;
  
  fwrite( &Nffa , sizeof( int ) , 1 , fw ) ;
  for( i = 0 ; i < Nffa ; i++ ) {
    fwrite( &Grid -> nstpx , sizeof( int ) , 1 , fw ) ;
    for( j = 0 ; j < Grid -> nstpx ; j++ ) {
      int nx = Grid -> nfx[j] ;
      fwrite( &Grid -> nstpy , sizeof( int ) , 1 , fw ) ;
      for( k = 0 ; k < Grid -> nstpy ; k++ ) {
	fwrite( &nx , sizeof( int ) , 1 , fw ) ;
	float buf[ nx ] ;
	int l ;
	for( l = 0 ; l < nx ; l++ ) {
	  buf[ l ] = (float)Grid -> Ffm[i][j][k][l] ;
	}
	fwrite( buf , sizeof( float ) , nx , fw ) ;
	DML_checksum_accum_crc32c( &cksumFfm_a , &cksumFfm_b , 
				   rank , (char*)buf , nx*sizeof(float) ) ;
	for( l = 0 ; l < nx ; l++ ) {
	  buf[ l ] = (float)Grid -> Ffp[i][j][k][l] ;
	}
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
  for( i = 0 ; i < (size_t)Grid -> NtayY ; i++ ) {
    fwrite( &Grid -> NY_tay , sizeof( int ) , 1 , fw ) ;
    fwrite( Grid -> TX[i] , sizeof( double ) , Grid -> NY_tay , fw  ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , Grid -> TX[i] ,
			       Grid->NY_tay * sizeof(double) ) ;
  }
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
  for( i = 0 ; i < (size_t)Grid -> NtayX ; i++ ) {
    fwrite( &Grid -> NX_tay , sizeof( int ) , 1 , fw ) ;
    fwrite( Grid -> TY[i] , sizeof( double ) , Grid -> NX_tay , fw  ) ;
    DML_checksum_accum_crc32c( &cksuma , &cksumb , 
			       i , Grid -> TY[i] ,
			       Grid->NX_tay * sizeof(double) ) ;
  }
  
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
