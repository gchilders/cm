/*

cm-impl.h - header file for internal use of the cm library

Copyright (C) 2009, 2010, 2012, 2015, 2016, 2018, 2021, 2022 Andreas Enge

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

#ifndef __CM_IMPL_H
#define __CM_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifdef WITH_MPI
#include <unistd.h> /* for usleep */
#include <mpi.h>
#endif
#include <pari/pari.h>
#include "cm.h"

#define CM_CLASS_DATADIR "."
#define CM_CLASS_TMPDIR "."

#define CM_FIELD_REAL    1
#define CM_FIELD_COMPLEX 2

#define CM_MIN(a,b) ((a) < (b) ? (a) : (b))
#define CM_MAX(a,b) ((a) > (b) ? (a) : (b))

typedef struct {
   cm_modular_t m;
   int_cl_t d;
      /* discriminant */
   int h12;
      /* number of forms up to negation */
   cm_form_t *form;
      /* primitive quadratic forms of discriminant d up to negation */
   ftype root;
      /* sqrt (-d); */
   ctype *eta;
      /* eta [i] contains the value of the Dedekind eta function in the
         root of the quadratic form form [i]. */
} cm_modclass_t;


typedef struct {
   unsigned long int counter [4];
   cm_timer_t timer [8];
}
__cm_stat_struct;

typedef __cm_stat_struct cm_stat_t [1];
typedef __cm_stat_struct *cm_stat_ptr;


#if defined (__cplusplus)
extern "C" {
#endif

/* version */
extern const char *cm_get_version (void);

/* number theoretic functions */
extern long int cm_nt_gcd (long int a, long int b);
extern int cm_nt_kronecker (int_cl_t a, int_cl_t b);
extern int cm_nt_is_prime (mpz_srcptr n);
extern unsigned long int cm_nt_next_prime (const unsigned long int n);
extern void cm_nt_factor (uint_cl_t d, uint_cl_t *factors,
   unsigned int *exponents);
extern uint_cl_t cm_nt_largest_factor (uint_cl_t n);
extern unsigned int cm_nt_mpz_tonelli_generator (mpz_ptr q, mpz_ptr z,
   mpz_srcptr p);
extern void cm_nt_mpz_tonelli_with_generator (mpz_ptr root,
   mpz_srcptr a, mpz_srcptr p, unsigned int e, mpz_srcptr q,
   mpz_srcptr z);
extern void cm_nt_mpz_tonelli_si_with_generator (mpz_ptr root,
   const long int a, mpz_srcptr p, unsigned int e, mpz_srcptr q,
   mpz_srcptr z);
extern void cm_nt_mpz_tonelli (mpz_ptr root, mpz_srcptr a, mpz_srcptr p);
extern void cm_nt_mpz_tonelli_si (mpz_ptr root, const long int a,
   mpz_srcptr p);
extern void cm_nt_mpz_tree_gcd (mpz_t *gcd, mpz_srcptr n, mpz_t *m,
   int no_m);
extern bool cm_nt_cget_zz (mpz_ptr out1, mpz_ptr out2, ctype in, ctype omega);

/* functions for computing q expansions of modular functions and addition
   chains */
extern void cm_qdev_init (cm_qdev_t *f, fprec_t prec);
extern void cm_qdev_clear (cm_qdev_t *f);
extern void cm_qdev_eval (ctype rop, cm_qdev_t f, ctype q1);
extern void cm_qdev_eval_fr (ftype rop, cm_qdev_t f, ftype q1);

/* function for evaluating modular functions */
extern void cm_modular_eta_series_fr (cm_modular_t m, ftype rop, ftype q_24);

/* functions depending on PARI */
extern void mpzx_oneroot_split_mod (mpz_ptr root, mpzx_srcptr f,
   mpz_srcptr p, bool verbose);
extern void cm_pari_oneroot (mpz_ptr root, mpzx_srcptr f, mpz_srcptr p);
extern mpz_t* cm_pari_find_roots (int *no, mpzx_srcptr f, mpz_srcptr p);
extern int cm_pari_classgroup (int_cl_t d, int_cl_t *ord, cm_form_t *gen);
extern int cm_pari_classgroup_2quotient (int_cl_t d, const int *p,
   int_cl_t *ord, cm_form_t *gen);

/* functions for integral polynomials */
extern void mpzx_init (mpzx_ptr f, int deg);
extern void mpzx_clear (mpzx_ptr f);
extern void mpzx_set_deg (mpzx_ptr f, int deg);
extern void mpzx_mod (mpzx_ptr g, mpzx_srcptr f, mpz_srcptr p);
extern bool cm_mpfrx_get_mpzx (mpzx_ptr g, mpfrx_srcptr f);
extern bool cm_mpcx_get_quadraticx (mpzx_ptr g, mpzx_ptr h, mpcx_srcptr f,
   int_cl_t d);
extern size_t mpzx_out_str (FILE* stream, int base, mpzx_srcptr f);
extern void mpzx_print_pari (FILE* file, mpzx_srcptr f, char *x);
extern void mpzxx_print_pari (FILE* file, mpzx_srcptr g, mpzx_srcptr h,
   char *x);

/* functions for number field towers */
extern void mpzx_tower_init (mpzx_tower_ptr twr, int levels, int *d);
extern void mpzx_tower_clear (mpzx_tower_ptr twr);
extern bool cm_mpfrx_tower_get_mpzx_tower (mpzx_tower_ptr tz,
   mpfrx_tower_srcptr tf);
extern bool cm_mpcx_tower_get_quadratic_tower (mpzx_tower_ptr t1,
   mpzx_tower_ptr t2, mpcx_tower_srcptr tc, int_cl_t d);
extern void mpzx_tower_print_pari (FILE* file, mpzx_tower_srcptr twr,
   char *fun, char *var);
extern void mpzxx_tower_print_pari (FILE* file, mpzx_tower_srcptr g,
   mpzx_tower_srcptr h, char *fun, char *var);


/* functions for classgroups of imaginary-quadratic number fields */

extern void cm_classgroup_init (cm_classgroup_t *cl, int_cl_t disc,
   const int *ramified, bool verbose);
extern void cm_classgroup_clear (cm_classgroup_t *cl);

extern void cm_classgroup_mpz_set_icl (mpz_t rop, int_cl_t op);
extern int_cl_t cm_classgroup_mpz_get_icl (mpz_t op);
extern uint_cl_t cm_classgroup_mod (int_cl_t a, uint_cl_t p);
extern int_cl_t cm_classgroup_gcd (int_cl_t a, int_cl_t b);

extern int_cl_t cm_classgroup_fundamental_primes (int_cl_t *primes,
      int_cl_t d);
extern int_cl_t cm_classgroup_fundamental_discriminant (int_cl_t d);

extern int_cl_t cm_classgroup_compute_c (int_cl_t a, int_cl_t b, int_cl_t d);
extern void cm_classgroup_reduce (cm_form_t *Q, int_cl_t d);
extern void cm_classgroup_compose (cm_form_t *Q, cm_form_t Q1,
   cm_form_t Q2, int_cl_t d);
extern void cm_classgroup_pow (cm_form_t *Q, cm_form_t P, uint_cl_t n,
   int_cl_t d);


/* functions for evaluating modular functions at quadratic integers via
   precomputations */

extern void cm_modclass_init (cm_modclass_t *mc, int_cl_t d, fprec_t prec,
   bool verbose);
extern void cm_modclass_clear (cm_modclass_t *mc);

extern void cm_modclass_eta_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_f_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, int e);
extern void cm_modclass_f1_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, int e);
extern void cm_modclass_gamma2_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_gamma3_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_j_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_multieta_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, const int *p, int e);
extern void cm_modclass_atkinhecke_level_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, unsigned long int l);


/* functions for class polynomials */
extern bool cm_class_write (cm_class_srcptr c, cm_param_srcptr param);
extern bool cm_class_read (cm_class_ptr c, cm_param_srcptr param);
extern mpz_t* cm_class_get_j_mod_p (int *no, cm_param_srcptr param,
   cm_class_srcptr c, mpz_srcptr p, const char* modpoldir, bool verbose);

/* functions for computing parameters of a complex multiplication curve */
extern void cm_curve_and_point_stat (mpz_ptr a, mpz_ptr b, mpz_ptr x,
   mpz_ptr y, cm_param_srcptr param, cm_class_srcptr c,
   mpz_srcptr p, mpz_srcptr l, mpz_srcptr co,
   const char* modpoldir, bool print, bool verbose, bool debug,
   cm_stat_t stat);


/* functions for ECPP */
extern void cm_pari_prime_product (mpz_ptr prim, unsigned long int a,
   unsigned long int b);
extern bool cm_pari_cornacchia (mpz_ptr t, mpz_ptr v, mpz_srcptr p,
   mpz_srcptr root, const int_cl_t d);
extern void cm_ecpp_compute_h_chunk (unsigned int *h, uint_cl_t Dmin,
   uint_cl_t Dmax);
extern mpz_t* cm_ecpp_compute_cardinalities (int *no_card,
   int_cl_t **card_d, int_cl_t *d, int no_d, mpz_srcptr N,
   long int *qstar, int no_qstar, mpz_t *qroot);
extern void cm_ecpp_one_step2 (mpz_t *cert2, mpz_t *cert1, int i,
   const char* modpoldir, bool verbose, bool debug, cm_stat_t stat);
extern bool cm_pari_ecpp_check (mpz_t **cert, int depth);


/* internal functions related to timers and counters */
extern void cm_stat_init (cm_stat_t stat);


/* functions operating on files */
extern bool cm_file_open_read_write (FILE **f, const char *filename);
extern bool cm_file_write_h (const char *tmpdir, const unsigned int *h,
   unsigned int n);
extern bool cm_file_read_h (const char *tmpdir, unsigned int *h,
   unsigned int n);
extern bool cm_file_write_primorial (const char *tmpdir, mpz_srcptr prim,
   const int i);
extern bool cm_file_read_primorial (const char *tmpdir, mpz_ptr prim,
   const int i);
extern bool cm_write_ecpp_cert1_line (FILE *f, mpz_t *line, cm_stat_t stat);
extern bool cm_write_ecpp_cert2_line (FILE *f, mpz_t *line, int no,
   cm_stat_t stat);
extern mpz_t** cm_file_read_ecpp_cert1 (int *depth, mpz_srcptr p,
   FILE *f, bool debug, cm_stat_t stat);
extern int cm_file_read_ecpp_cert2 (mpz_t **c, mpz_srcptr p,
   FILE *f, bool debug, cm_stat_t stat);
extern void cm_file_write_ecpp_cert_pari (FILE *f, mpz_t **c, int l);
extern void cm_file_write_ecpp_cert_primo (FILE *f, mpz_t **c, int l);

/* functions for MPI */
extern void cm_mpi_queue_push (int rank);
extern int cm_mpi_queue_pop (void);
extern unsigned long int cm_mpi_compute_B (void);
extern void cm_mpi_broadcast_init (bool verbose, bool debug);
extern void cm_mpi_broadcast_N (mpz_srcptr N);
extern void cm_mpi_broadcast_sqrt (int no_qstar, long int *qstar,
   mpz_t *qroot);
extern void cm_mpi_clear_N (void);
extern void cm_mpi_submit_primorial (char *tmpdir);
extern void cm_mpi_get_primorial (int rank, double *t);
extern void cm_mpi_submit_tonelli (int rank, int job, const long int a);
extern void cm_mpi_get_tonelli (mpz_ptr root, int rank, double *t);
extern void cm_mpi_submit_ecpp_one_step2 (int rank, int job, mpz_t *cert1,
   const char* modpoldir);
extern void cm_mpi_get_ecpp_one_step2 (mpz_t *cert2, int rank,
   cm_stat_ptr stat);
extern void cm_mpi_submit_curve_cardinalities (int rank, int job,
   int_cl_t *d, int no_d);
extern mpz_t* cm_mpi_get_curve_cardinalities (int *no_card,
   int_cl_t **card_d, int rank, double *t);
extern void cm_mpi_submit_is_prime (int rank, int job, mpz_srcptr n);
extern bool cm_mpi_get_is_prime (int rank, double *t);
extern void cm_mpi_submit_h_chunk (int rank, int job, uint_cl_t Dmin,
   uint_cl_t Dmax);
extern void cm_mpi_get_h_chunk (unsigned int *h, int rank, double *t);
extern void cm_mpi_submit_tree_gcd (mpz_t *n, int no_n);
extern void cm_mpi_get_tree_gcd (mpz_t *gcd, int no_n, double *t);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_IMPL_H */
