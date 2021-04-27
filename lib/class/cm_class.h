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
   int h;
      /* the class number */
   int minpoly_deg;
      /* the degree of the minimal polynomial; usually h */
   mpz_t *minpoly;
      /* real part of the minimal polynomial of the function over Q             */
   mpz_t *minpoly_complex;
      /* Only meaningful in the complex case; then the minimal polynomial is    */
      /* decomposed into two parts over the integral basis                      */
      /* [1, sqrt (D)/2] resp. [1, (1 + sqrt (D))/2]; the first part is in      */
      /* minpoly, the second one in this variable.                              */
   int tower_levels;
   int *tower_d;
   mpz_t ***tower;
      /* These fields are meaningful only when the class field is
         decomposed as a tower; they represent the number of levels in the
         tower, the degree sequence (from bottom to top) and the
         polynomials defining the extension. Hereby, tower [i][j] encodes
         the rounded polynomial W [i][j] in an mpfrx_tower_t, with
         tower [i][j][k] being the coefficient of degree k. The degrees are
         not stored separately, but can be derived from tower_d: The degree
         of tower [0][0] is d [0] (with leading coefficient 1 that is also
         stored), that of tower [i][j] for i >= 1 is
         d [0] * ... * d [i-1] - 1, with potentially leading coefficients 0
         that are also stored. */
} cm_class_t;


#if defined (__cplusplus)
extern "C" {
#endif

/* functions for class polynomials */
extern void cm_class_init (cm_class_t *c, int_cl_t d, char inv,
   bool verbose);
extern void cm_class_clear (cm_class_t *c);
extern void cm_class_compute_minpoly (cm_class_t c, bool tower,
   bool checkpoints, bool disk, bool print, bool verbose);

/* functions for computing parameters of a complex multiplication curve      */
extern void cm_curve_compute_curve (int_cl_t d, char inv, int fieldsize,
   const char* modpoldir, bool readwrite, bool print, bool verbose);
#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_CLASS_H */
