/*

cm_class-impl.h - header file for internal use of the cm_class library

Copyright (C) 2009, 2010, 2015, 2016 Andreas Enge

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

#ifndef __CM_CLASS_IMPL_H
#define __CM_CLASS_IMPL_H

#include "mpfpx.h"
#include "cm_class.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define CM_CLASS_DATADIR "."
#define CM_CLASS_TMPDIR "."

#define CM_FIELD_REAL    1
#define CM_FIELD_COMPLEX 2


typedef enum {real, complex, drop}
   cm_embedding_t;
   /* When applied to class group entries, "real" indicates an ambiguous     */
   /* form, "complex" a non-ambiguous form with b > 0 and "conj" a           */
   /* non-ambiguous form with b < 0. The terminology is explained by the     */
   /* fact that for j, these forms yield real and complex values, with       */
   /* "conj" being the complex conjugate of a corresponding "complex".       */
   /* When applied to an N-system for other functions than j,                */
   /* "real" and "complex" indicate a real or complex conjugate. "drop" is   */
   /* used only in the case of a real class polynomial; it means that the    */
   /* conjugate is the complex conjugate of another one and can thus be      */
   /* dropped.                                                               */

typedef struct {
   int_cl_t a, b;
   cm_embedding_t emb;
} cm_form_t;

typedef struct cm_avl_t {
   struct cm_avl_t *l, *r; /* left and right subtrees */
   signed char     b; /* balance factor; -1, 0 or 1 */
   cm_form_t       c; /* content */
} cm_avl_t;

typedef struct {
   int_cl_t d;
   cm_form_t *form;
      /* contains a set of representatives of quadratic forms of             */
      /* discriminant d. Only forms with b >= 0 are stored, since the others */
      /* simply yield complex conjugate eta values.                          */
   int h1, h2, h12, h;
      /* the numbers of ambiguous reduced forms, pairs of non-ambiguous      */
      /* reduced forms, h1+h2 and the class number                           */
} cm_classgroup_t;

typedef struct {
   cm_modular_t m;
   cm_classgroup_t cl;
   ftype root;
      /* sqrt (-cl.d); */
   ftype sqrt2_over2, sqrt2_over4;
   ctype *eta;
      /* contains the values of eta with respect to the entries of           */
      /* cl.form with the same row index. So only cl.h12 values are stored.  */
   cm_classgroup_t cl2;
   ctype *eta2;
   ftype root2;
      /* Space for a second class group and associated eta values; intended  */
      /* to be filled with an order of smaller conductor. Currently, is is   */
      /* used for the Weber function with odd discriminant D; then cl.d==4*D */
      /* and cl2.d==D.                                                       */
} cm_modclass_t;


#if defined (__cplusplus)
extern "C" {
#endif

/* functions depending on PARI                                               */
extern void cm_pari_onefactor (mpz_t *res, mpz_t *f, int deg, int deg_factor,
   mpz_t p, bool verbose);
extern void cm_pari_oneroot (mpz_t root, mpz_t *f, int deg, mpz_t p,
   bool verbose);
extern mpz_t* cm_pari_find_roots (mpz_t *f, int deg, mpz_t p, int *no);


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
extern int cm_classgroup_h (int *h1, int *h2, int_cl_t d);

extern bool cm_classgroup_avl_insert (cm_avl_t **t, cm_form_t c);
extern int cm_classgroup_avl_count (cm_avl_t *t);
extern void cm_classgroup_avl_flatten (cm_form_t **list, cm_avl_t *t);
extern void cm_classgroup_avl_delete (cm_avl_t *t);

extern int_cl_t cm_classgroup_compute_c (int_cl_t a, int_cl_t b, int_cl_t d);
extern void cm_classgroup_reduce (cm_form_t *Q, int_cl_t d);
extern void cm_classgroup_compose (cm_form_t *Q, cm_form_t Q1,
   cm_form_t Q2, int_cl_t d);
extern cm_form_t cm_classgroup_prime_form (int_cl_t p, int_cl_t d);


/* functions for evaluating modular functions at quadratic integers via
   precomputations */

extern void cm_modclass_init (cm_modclass_t *mc, cm_classgroup_t cl,
   cm_classgroup_t cl2, fprec_t prec, bool checkpoints, bool verbose);
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
extern void cm_class_write (cm_class_t c);
extern bool cm_class_read (cm_class_t c);

extern mpz_t* cm_class_get_j_mod_P (int_cl_t d, char inv, mpz_t P, int *no,
   const char* modpoldir, bool readwrite, bool verbose);


#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_CLASS_IMPL_H */
