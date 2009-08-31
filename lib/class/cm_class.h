#ifndef __CM_CLASS_H
#define __CM_CLASS_H

#include "cm_common.h"
#include <inttypes.h>

#define CM_INVARIANT_NONE      '\0'
#define CM_INVARIANT_J         'j'
#define CM_INVARIANT_GAMMA2    '2'
#define CM_INVARIANT_GAMMA3    '3'
#define CM_INVARIANT_WEBER     'w'
#define CM_INVARIANT_DOUBLEETA 'd'
#define CM_INVARIANT_SIMPLEETA 's'
#define CM_INVARIANT_ATKIN     'a'


typedef int_fast64_t int_cl_t;
typedef uint_fast64_t uint_cl_t;
   /* integer types used to represent discriminants, conductors and class    */
   /* group entries  */
#define PRIicl PRIiFAST64
#define PRIucl PRIuFAST64
#define SCNicl SCNiFAST64


#if defined (__cplusplus)
extern "C" {
#endif

/* functions for computing parameters of a complex multiplication curve      */
extern void cm_curve_compute_curve (int_cl_t d, char inv,
   const char* modpoldir, bool verbose);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_CLASS_H */
