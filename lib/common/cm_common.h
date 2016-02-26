/*

cm_common.h - header file for the cm_common library

Copyright (C) 2009, 2010, 2012, 2015, 2016 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifndef __CM_COMMON_H
#define __CM_COMMON_H

#include <sys/times.h>
#include <stdbool.h>
#include <zlib.h>
#include "cm_arith.h"
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
extern int cm_nt_is_prime_l (const unsigned long int prime);
extern unsigned long int cm_nt_next_prime (const unsigned long int n);
extern long int cm_nt_gcd (long int a, long int b);
extern long int cm_nt_gcdext (long int *u, long int *v, long int a,
   long int b);
extern int cm_nt_kronecker (long int a, long int b);
extern long int cm_nt_sqrt (const unsigned long int n);
extern void cm_nt_factor (long int d, unsigned long int *factors,
   unsigned int *exponents);

extern void cm_nt_mpz_tonelli_z (mpz_t root, mpz_t a, mpz_t p);
extern void cm_nt_mpz_tonelli (mpz_t root, const long int a, mpz_t p);

extern void cm_nt_elliptic_curve_multiply (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t m, mpz_t a, mpz_t p);
extern void cm_nt_elliptic_curve_random (mpz_t P_x, mpz_t P_y,
   mpz_t cofactor, mpz_t a, mpz_t b, mpz_t p);

extern bool cm_nt_fget_z (mpz_t out, ftype in);

/* functions for evaluating modular functions */
extern void cm_modular_fundamental_domain (cptr z);
extern void cm_modular_init (cm_modular_t *m, fprec_t prec);
extern void cm_modular_clear (cm_modular_t *m);
extern void cm_modular_eta_transform (cm_modular_t m, ctype rop, ctype z,
   cm_matrix_t M);
extern void cm_modular_eta_series (cm_modular_t m, ctype rop, ctype q_24);
extern void cm_modular_eta_eval (cm_modular_t m, ctype rop, ctype op);
extern void cm_modular_eta_eval_fr (cm_modular_t m, ftype rop, ftype op);
extern void cm_modular_atkinhecke_eval (cm_modular_t m, ctype rop, ctype op,
   unsigned long int l, unsigned long int r);
extern void cm_modular_atkinhecke_level_eval (cm_modular_t m, ctype rop,
   ctype op, unsigned long int l);

/* functions for evaluating modular functions using the AGM */
extern void cm_fem_eta_eval (cm_modular_t m, ctype rop, ctype op);

/* functions reading modular polynomials */
extern mpz_t* cm_modpol_read_specialised_mod (int* n, int level, char type,
   mpz_t p, mpz_t x, const char * datadir);
extern void cm_modpol_print_pari (int level, char type, const char* datadir);
extern void cm_modpol_print_magma (int level, char type, const char* datadir);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_COMMON_H */
