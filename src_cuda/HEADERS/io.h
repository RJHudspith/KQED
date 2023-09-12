#ifndef IO_H
#define IO_H

int
read_ff( struct Grid_coeffs *Grid ) ;

int
read_TAYLORX( struct Grid_coeffs *Grid ) ;

int
read_TAYLORY( struct Grid_coeffs *Grid ) ;

void
write_ff( const struct Grid_coeffs *Grid ) ;

void
write_TAYLORX( const struct Grid_coeffs *Grid ) ;

void
write_TAYLORY( const struct Grid_coeffs *Grid ) ;

#endif
