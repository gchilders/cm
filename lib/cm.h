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

#include <time.h>
#include <stdbool.h>
#include <zlib.h>
#include <inttypes.h>
#include <sys/time.h>
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

#define CM_SUBFIELD_NEVER     0
#define CM_SUBFIELD_PREFERRED 1
#define CM_SUBFIELD_OPTIMAL   2


typedef struct {
   clock_t time_old;
   double elapsed;
   struct timeval wc_time_old [1];
   double wc_elapsed;

} __cm_timer_struct;
typedef __cm_timer_struct cm_timer_t [1];

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


/* Types for quadratic forms and class groups */
typedef struct {
   int_cl_t a, b;
} cm_form_t;

typedef struct {
   int_cl_t d;
   cm_form_t *form;
      /* contains a set of representatives of quadratic forms of             */
      /* discriminant d.                                                     */
   int h;
      /* the class number */
   int levels;
   int *deg;
      /* Number of entries and their cardinalities (from back to front) in
         the normal series associated to the ordering of the quadratic
         forms, or equivalently the sequence of degrees (from bottom to
         top) in a Galois tower decomposition of the class field. */
} cm_classgroup_t;


/* Types for CM parameters. */

typedef struct {
   int_cl_t d;
      /* The (not necessarily fundamental) discriminant of an imaginary-
         quadratic order. */
   char invariant;
      /* Constant describing which type of invariant is used. */
   int field;
      /* A constant describing whether we are working over the real or the
         complex numbers. */
   int p [6], e, s;
      /* The actual parameters attached to the class invariant.
         p is a 0-terminated list of integers (often the primes dividing the
         level); s is the canonical power, e the power actually used. */
   int r [6];
      /* The 0-terminated list of ramified primes not dividing the
         conductor, which is used for double and multiple eta quotients
         to determine subfields of the ring class field. */
   char str [255];
      /* A string encoding the previous characters, used in files and their
         names. */
} __cm_param_struct;

typedef __cm_param_struct cm_param_t [1];
typedef __cm_param_struct *cm_param_ptr;
typedef const __cm_param_struct *cm_param_srcptr;


/* Type for class polynomials and class field towers. */

typedef struct {
   cm_classgroup_t cl;
      /* The class group; it also contains the discriminant d and the class
         number h. */
   int_cl_t dfund;
      /* The fundamental discriminant attached to d, needed for rounding to
         quadratic integers in the complex case. */
   mpzx_t classpol;
      /* Real part of the class polynomial of the function over Q. */
   mpzx_t classpol_c;
      /* Only meaningful in the complex case; then the minimal polynomial is
         decomposed into two parts over the integral basis
         [1, sqrt (D)/2] resp. [1, (1 + sqrt (D))/2]; the first part is in
         classpol, the second one in this field. */
   mpzx_tower_t tower;
      /* This field is meaningful only when the class field is decomposed
         as a tower; it represents the polynomials defining the extensions,
         in the same format as an mpfrx_tower or an mpcx_tower. */
   mpzx_tower_t tower_c;
      /* This field is meaningful only in the complex case and when the
         class field is decomposed as a tower; it contains the entries of
         the defining polynomials in the second element of the integral
         basis as explained for classpol_c. */
   int field;
      /* This is a duplicate of the field with the same name in cm_param,
         but it makes the structure self-contained with respect to which
         fields are initialised. */
   bool computed_classpol;
   bool computed_tower;
      /* These fields store whether the class polynomial or the tower
         decomposition are stored in the variable. */
} __cm_class_struct;

typedef __cm_class_struct cm_class_t [1];
typedef __cm_class_struct *cm_class_ptr;
typedef const __cm_class_struct *cm_class_srcptr;

#if defined (__cplusplus)
extern "C" {
#endif

/* functions for measuring the passing time */
extern void cm_timer_start (cm_timer_t t);
extern void cm_timer_reset (cm_timer_t t);
extern void cm_timer_continue (cm_timer_t t);
extern void cm_timer_stop (cm_timer_t t);
extern double cm_timer_get (cm_timer_t t);
extern double cm_timer_wc_get (cm_timer_t t);

/* generic functions for opening files */
extern bool cm_file_open_write (FILE **f, char *filename);
extern bool cm_file_open_read (FILE **f, char *filename);
extern void cm_file_close (FILE *f);
extern void cm_file_gzopen_write (gzFile *f, char *filename);
extern void cm_file_gzopen_read (gzFile *f, char *filename);
extern void cm_file_gzclose (gzFile f);

/* functions for rounding floating point numbers to rational or quadratic
   integers */
extern bool cm_nt_fget_z (mpz_t out, ftype in);

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
extern void cm_modpol_read_specialised_mod (mpzx_ptr pol, int level,
   char type, mpz_srcptr p, mpz_srcptr x, const char * datadir);
extern void cm_modpol_print_pari (int level, char type, const char* datadir);
extern void cm_modpol_print_magma (int level, char type, const char* datadir);

/* functions concerned with CM parameters */
extern bool cm_param_init (cm_param_ptr param, int_cl_t d, char invariant,
   int maxdeg, int subfield, bool verbose);

/* functions depending on PARI */
extern void cm_pari_init (void);
extern void cm_pari_clear (void);
extern bool cm_pari_eval_int (mpz_ptr n, char *e);

/* functions for class polynomials */
extern void cm_class_init (cm_class_ptr c, cm_param_srcptr param,
   bool verbose);
extern void cm_class_clear (cm_class_ptr c);
extern double cm_class_height_factor (cm_param_srcptr param);
extern bool cm_class_compute (cm_class_ptr c, cm_param_srcptr param,
   bool classpol, bool tower, bool verbose);
extern void cm_class_print_pari (FILE* file, cm_class_srcptr c,
   char *fun, char *var);

/* functions for computing parameters of a complex multiplication curve */
extern void cm_curve_crypto_param (mpz_ptr p, mpz_ptr n, mpz_ptr l,
   mpz_ptr c, int_cl_t d, int fieldsize, bool verbose);
extern void cm_curve_and_point (mpz_ptr a, mpz_ptr b, mpz_ptr x, mpz_ptr y,
   cm_param_srcptr param, cm_class_srcptr c,
   mpz_srcptr p, mpz_srcptr l, mpz_srcptr co,
   const char* modpoldir, bool print, bool verbose);

/* functions for ECPP */
extern bool cm_ecpp (mpz_srcptr N, const char* modpoldir, bool tower,
   bool print, bool check, bool verbose, bool debug);

/* functions for MPI */
void cm_mpi_init (void);
void cm_mpi_clear (void);

#if defined (__cplusplus)
}
#endif

#endif /* ifndef __CM_H */
