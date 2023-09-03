#ifndef INIT_H
#define INIT_H

__host__
void
free_QED_temps( struct QED_kernel_temps *t ) ;

__host__
int
initialise( struct QED_kernel_temps *t ) ;

#endif
