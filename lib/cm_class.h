/*

cm_class.h - header file for the cm_class library

Copyright (C) 2009, 2010, 2018, 2021 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

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
#define CM_INVARIANT_MULTIETA  'm'
#define CM_INVARIANT_ATKIN     'a'


typedef int_fast64_t int_cl_t;
typedef uint_fast64_t uint_cl_t;
   /* integer types used to represent discriminants, conductors and class    */
   /* group entries  */
#define PRIicl PRIiFAST64
#define PRIucl PRIuFAST64
#define SCNicl SCNiFAST64


/* Type for polynomials, modelled after mpfrx and mpcx, except that the
   size and the degree are not handled separately; this does not matter
   for our application, where all variables are essentially "final" and
   only assigned once by rounding from a floating point polynomial. */
typedef struct {
   int deg;
      /* A degree of -1 indicates the zero polynomial. */
   mpz_t *coeff;
}
__mpzx_struct;

typedef __mpzx_struct mpzx_t[1];
typedef __mpzx_struct *mpzx_ptr;
typedef const __mpzx_struct *mpzx_srcptr;


/* Type for the definition of number field towers, modelled after
   mpfrx_tower and mpcx_tower. */
typedef struct {
   int levels;
   int* d;
   int deg;
   mpzx_t** W;
}
__mpzx_tower_struct;

typedef __mpzx_tower_struct mpzx_tower_t [1];
typedef __mpzx_tower_struct *mpzx_tower_ptr;
typedef const __mpzx_tower_struct *mpzx_tower_srcptr;


typedef struct {
   char invariant;
      /* a constant describing which invariant is actually used                 */
   int field;
      /* a constant describing whether we are working over the real or the      */
      /* complex numbers                                                        */
   int p [6], e, s;
      /* some parameters of the class invariant                                 */
      /* p is a 0-terminated list of integers (often the primes dividing the    */
      /* level); s is the canonical power, e the power actually used.           */
   char paramstr [255];
      /* a string encoding the previous characters, used in files and their     */
      /* names                                                                  */
   int_cl_t d;
      /* the discriminant                                                       */
   int_cl_t dfund;
      /* The fundamental discriminant attached to d, needed for rounding to
         quadratic integers in the complex case. */
   int h;
      /* the class number */
   mpzx_t minpoly;
      /* real part of the minimal polynomial of the function over Q             */
   mpzx_t minpoly_complex;
      /* Only meaningful in the complex case; then the minimal polynomial is    */
      /* decomposed into two parts over the integral basis                      */
      /* [1, sqrt (D)/2] resp. [1, (1 + sqrt (D))/2]; the first part is in      */
      /* minpoly, the second one in this variable.                              */
   mpzx_tower_t tower;
      /* This field is meaningful only when the class field is decomposed
         as a tower; it represents the polynomials defining the extensions,
         in the same format as an mpfrx_tower or an mpcx_tower. */
   mpzx_tower_t tower_complex;
      /* This field is meaningful only in the complex case and when the
         class field is decomposed as a tower; it contains the entries of
         the defining polynomials in the second element of the integral
         basis as explained for minpoly_complex. */
} cm_class_t;


#if defined (__cplusplus)
extern "C" {
#endif

/* functions for class polynomials */
extern void cm_class_init (cm_class_t *c, int_cl_t d, char inv,
   bool verbose);
extern void cm_class_clear (cm_class_t *c);
extern bool cm_class_compute_minpoly (cm_class_t c, bool tower,
   bool checkpoints, bool disk, bool print, bool verbose);

/* functions for computing parameters of a complex multiplication curve      */
extern void cm_curve_compute_curve (int_cl_t d, char inv, int fieldsize,
   const char* modpoldir, bool readwrite, bool print, bool verbose);
#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_CLASS_H */