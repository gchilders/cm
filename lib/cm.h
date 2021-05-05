/*

cm.h - header file for the cm library

Copyright (C) 2009, 2010, 2012, 2015, 2016, 2018, 2021 Andreas Enge

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

#ifndef __CM_H
#define __CM_H

#include <sys/times.h>
#include <stdbool.h>
#include <zlib.h>
#include <inttypes.h>
#include "cm-arith.h"
#include "mpfrcx.h"


#define CM_MODPOL_J           'j'
#define CM_MODPOL_DOUBLEETA   'd'
#define CM_MODPOL_SIMPLEETA   's'
#define CM_MODPOL_MULTIETA    'm'
#define CM_MODPOL_WEBER       'w'
#define CM_MODPOL_WEBERSQUARE '2'
#define CM_MODPOL_U           'u'
#define CM_MODPOL_ATKIN       'a'
#define CM_MODPOL_W35         '5'
#define CM_MODPOL_W39         '3'
#define CM_MODPOL_GAMMA2      'g'
#define CM_MODPOL_ATKIN47     '4'
#define CM_MODPOL_ATKIN59     'A'
#define CM_MODPOL_ATKIN71     '7'
#define CM_MODPOL_RAMANUJAN   'r'
#define CM_MODPOL_MULTI30     '0'
#define CM_MODPOL_MULTI42     'x'
#define CM_MODPOL_MULTI70     'y'
#define CM_MODPOL_MULTI78     '8'

#define CM_INVARIANT_NONE      '\0'
#define CM_INVARIANT_J         'j'
#define CM_INVARIANT_GAMMA2    '2'
#define CM_INVARIANT_GAMMA3    '3'
#define CM_INVARIANT_WEBER     'w'
#define CM_INVARIANT_DOUBLEETA 'd'
#define CM_INVARIANT_SIMPLEETA 's'
#define CM_INVARIANT_MULTIETA  'm'
#define CM_INVARIANT_ATKIN     'a'


typedef struct {
   struct tms time_old;
   double     elapsed;
} __cm_timer_struct;
typedef __cm_timer_struct cm_timer [1];

typedef struct {
   long int a, b, c, d;
} cm_matrix_t;

typedef struct {
   long int **chain;
      /* data structure for holding addition chains                          */
      /* entry 0: the value of the exponent; chain [0][0] must be 0 and      */
      /*                                     chain [1][0] must be 1          */
      /* entry 1: the rule for obtaining this exponent, with the following   */
      /*          meaning:                                                   */
      /*          1: 2*i1                                                    */
      /*          2: i1 + i2                                                 */
      /*          3: 2*i1 + i2                                             */
      /* entries 2 to 3: the indices i1 and i2 yielding this exponent        */
      /* entry 4: the coefficient with which the term contributes to the     */
      /*          function (0 if it is only used as an auxiliary term)       */
   int length;
      /* the number of terms actually computed for the addition chain */
} cm_qdev_t;

typedef struct {
   fprec_t prec;
   ctype zeta48inv;
   ftype pi;
   ctype log_zeta24;
   ctype twopii;
   ctype zeta24 [24];
   ftype sqrt2;
   cm_qdev_t eta;
} cm_modular_t;


/* Integer types used to represent discriminants, conductors and class
   group entries  */
typedef int_fast64_t int_cl_t;
typedef uint_fast64_t uint_cl_t;
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

/* functions for measuring the passing time */
extern void cm_timer_start (cm_timer clock);
extern void cm_timer_stop (cm_timer clock);
extern double cm_timer_get (cm_timer clock);

/* generic functions for opening files */
extern bool cm_file_open_write (FILE **f, char *filename);
extern bool cm_file_open_read (FILE **f, char *filename);
extern void cm_file_close (FILE *f);
extern void cm_file_gzopen_write (gzFile *f, char *filename);
extern void cm_file_gzopen_read (gzFile *f, char *filename);
extern void cm_file_gzclose (gzFile f);

/* different functions for number theoretic computations */
extern int cm_nt_is_prime (mpz_t a);
extern unsigned long int cm_nt_next_prime (const unsigned long int n);
extern long int cm_nt_gcd (long int a, long int b);
extern int cm_nt_kronecker (long int a, long int b);
extern void cm_nt_factor (long int d, unsigned long int *factors,
   unsigned int *exponents);

extern void cm_nt_mpz_tonelli_z (mpz_t root, mpz_t a, mpz_t p);
extern void cm_nt_mpz_tonelli (mpz_t root, const long int a, mpz_t p);

extern void cm_nt_elliptic_curve_multiply (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t m, mpz_t a, mpz_t p);
extern void cm_nt_elliptic_curve_random (mpz_t P_x, mpz_t P_y,
   mpz_t cofactor, mpz_t a, mpz_t b, mpz_t p);

/* functions for rounding floating point numbers to rational or quadratic
   integers */
extern bool cm_nt_fget_z (mpz_t out, ftype in);
extern bool cm_nt_cget_zz (mpz_ptr out1, mpz_ptr out2, ctype in, ctype omega);

/* functions for evaluating modular functions */
extern void cm_modular_fundamental_domain (cptr z);
extern void cm_modular_init (cm_modular_t *m, fprec_t prec);
extern void cm_modular_clear (cm_modular_t *m);
extern int cm_modular_eta_transform (long int *e, ctype czplusd, ctype z,
   cm_matrix_t M);
extern void cm_modular_eta_series (cm_modular_t m, ctype rop, ctype q_24);
extern void cm_modular_eta_eval (cm_modular_t m, ctype rop, ctype op);
extern void cm_modular_eta_eval_fr (cm_modular_t m, ftype rop, ftype op);
extern void cm_modular_atkinhecke_eval (cm_modular_t m, ctype rop, ctype op,
   unsigned long int l, unsigned long int r);
extern void cm_modular_atkinhecke_level_eval (cm_modular_t m, ctype rop,
   ctype op, unsigned long int l);

/* functions reading modular polynomials */
extern mpz_t* cm_modpol_read_specialised_mod (int* n, int level, char type,
   mpz_t p, mpz_t x, const char * datadir);
extern void cm_modpol_print_pari (int level, char type, const char* datadir);
extern void cm_modpol_print_magma (int level, char type, const char* datadir);

/* functions for class polynomials */
extern void cm_class_init (cm_class_t *c, int_cl_t d, char inv,
   bool verbose);
extern void cm_class_clear (cm_class_t *c);
extern bool cm_class_compute_minpoly (cm_class_t c, bool tower,
   bool disk, bool print, bool verbose);

/* functions for computing parameters of a complex multiplication curve      */
extern void cm_curve_compute_curve (int_cl_t d, char inv, int fieldsize,
   const char* modpoldir, bool readwrite, bool print, bool verbose);

#if defined (__cplusplus)
}
#endif

#endif /* ifndef __CM_H */
