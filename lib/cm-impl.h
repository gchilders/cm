/*

cm-impl.h - header file for internal use of the cm library

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

#ifndef __CM_IMPL_H
#define __CM_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <pari/pari.h>
#include "cm.h"

#define CM_CLASS_DATADIR "."
#define CM_CLASS_TMPDIR "."

#define CM_FIELD_REAL    1
#define CM_FIELD_COMPLEX 2


typedef struct {
   int_cl_t a, b;
} cm_form_t;


typedef struct {
   int_cl_t d;
   cm_form_t *form;
      /* contains a set of representatives of quadratic forms of             */
      /* discriminant d.                                                     */
   int *conj;
      /* references how the forms in form are related by "conjugation",
         that is, negation; conj [i] == j if and only if form [j] is the
         inverse of form [i] */
   int h;
      /* the class number */
   int levels;
   int *deg;
      /* Number of entries and their cardinalities (from back to front) in
         the normal series associated to the ordering of the quadratic
         forms, or equivalently the sequence of degrees (from bottom to
         top) in a Galois tower decomposition of the class field. */
} cm_classgroup_t;


typedef struct {
   cm_modular_t m;
   cm_classgroup_t cl;
   ftype root;
      /* sqrt (-cl.d); */
   ftype sqrt2_over2, sqrt2_over4;
   ctype *eta;
      /* contains the values of eta with respect to the entries of
         cl.form with the same row index. */
} cm_modclass_t;


#if defined (__cplusplus)
extern "C" {
#endif

/* number theoretic functions */
extern void cm_nt_mpz_tonelli_z (mpz_ptr root, mpz_srcptr a, mpz_srcptr p);
extern void cm_nt_mpz_tonelli (mpz_ptr root, const long int a, mpz_srcptr p);
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
extern void cm_pari_oneroot (mpz_ptr root, mpzx_srcptr f, mpz_srcptr p,
   bool verbose);
extern mpz_t* cm_pari_find_roots (int *no, mpzx_srcptr f, mpz_srcptr p);
extern int cm_pari_classgroup (int_cl_t d, int_cl_t *ord, cm_form_t *gen);

/* functions for integral polynomials */
extern void mpzx_init (mpzx_ptr f, int deg);
extern void mpzx_clear (mpzx_ptr f);
extern bool cm_mpfrx_get_mpzx (mpzx_ptr g, mpfrx_srcptr f);
extern bool cm_mpcx_get_quadraticx (mpzx_ptr g, mpzx_ptr h, mpcx_srcptr f,
   int_cl_t d);
extern size_t mpzx_out_str (FILE* stream, int base, mpzx_srcptr f);
extern void mpzx_print_pari (FILE* file, mpzx_srcptr f, char *x);

/* functions for number field towers */
extern void mpzx_tower_init (mpzx_tower_ptr twr, int levels, int *d);
extern void mpzx_tower_clear (mpzx_tower_ptr twr);
extern bool cm_mpfrx_tower_get_mpzx_tower (mpzx_tower_ptr tz,
   mpfrx_tower_srcptr tf);
extern bool cm_mpcx_tower_get_quadratic_tower (mpzx_tower_ptr t1,
   mpzx_tower_ptr t2, mpcx_tower_srcptr tc, int_cl_t d);
extern void mpzx_tower_print_pari (FILE* file, mpzx_tower_srcptr twr,
   char *fun, char *var);


/* functions for classgroups of imaginary-quadratic number fields */

extern void cm_classgroup_init (cm_classgroup_t *cl, int_cl_t disc,
   bool verbose);
extern void cm_classgroup_clear (cm_classgroup_t *cl);

extern void cm_classgroup_mpz_set_icl (mpz_t rop, int_cl_t op);
extern int_cl_t cm_classgroup_mpz_get_icl (mpz_t op);
extern uint_cl_t cm_classgroup_mod (int_cl_t a, uint_cl_t p);
extern int_cl_t cm_classgroup_gcd (int_cl_t a, int_cl_t b);
extern int cm_classgroup_kronecker (int_cl_t a, int_cl_t b);

extern void cm_classgroup_factor (int_cl_t d,
      uint_cl_t *factors, unsigned int *exponents);
extern int_cl_t cm_classgroup_fundamental_discriminant (int_cl_t d);

extern int_cl_t cm_classgroup_compute_c (int_cl_t a, int_cl_t b, int_cl_t d);
extern void cm_classgroup_reduce (cm_form_t *Q, int_cl_t d);
extern void cm_classgroup_compose (cm_form_t *Q, cm_form_t Q1,
   cm_form_t Q2, int_cl_t d);
extern void cm_classgroup_pow (cm_form_t *Q, cm_form_t P, uint_cl_t n,
   int_cl_t d);

extern int cm_classgroup_normalseries (int_cl_t disc, int_cl_t *ord,
   cm_form_t *gen);


/* functions for evaluating modular functions at quadratic integers via
   precomputations */

extern void cm_modclass_init (cm_modclass_t *mc, cm_classgroup_t cl,
   fprec_t prec, bool verbose);
extern void cm_modclass_clear (cm_modclass_t *mc);

extern void cm_modclass_eta_eval_quad (ctype rop, cm_modular_t m,
   cm_classgroup_t cl, ctype *eta, int_cl_t a, int_cl_t b, ftype root);
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
   int_cl_t a, int_cl_t b, int *p, int e);
extern void cm_modclass_atkinhecke_level_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, unsigned long int l);


/* functions for class polynomials */
extern bool cm_class_write (cm_class_t c);
extern bool cm_class_read (cm_class_t c);

extern mpz_t* cm_class_get_j_mod_P (int_cl_t d, char inv, mpz_t P, int *no,
   const char* modpoldir, bool pari, bool readwrite, bool tower,
   bool verbose);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_IMPL_H */
