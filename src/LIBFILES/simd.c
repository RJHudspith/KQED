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
}
