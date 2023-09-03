#ifndef IO_H
#define IO_H

__host__
int
read_ff( struct Grid_coeffs *Grid ) ;

__host__
int
read_TAYLORX( struct Grid_coeffs *Grid ) ;

__host__
int
read_TAYLORY( struct Grid_coeffs *Grid ) ;

__host__
void
write_ff( const struct Grid_coeffs *Grid ) ;

__host__
void
write_TAYLORX( const struct Grid_coeffs *Grid ) ;

__host__
void
write_TAYLORY( const struct Grid_coeffs *Grid ) ;

#endif
