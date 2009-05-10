#ifndef __CM_CLASS_IMPL_H
#define __CM_CLASS_IMPL_H

#include "mpfpx.h"
#include "cm_class.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define CM_CLASS_DATADIR "."
#define CM_CLASS_TMPDIR "."

#define CM_FIELD_REAL    1
#define CM_FIELD_COMPLEX 2


typedef struct {
   int_cl_t d;
   int_cl_t **form;
      /* contains a set of representatives of quadratic forms of             */
      /* discriminant d. Each row corresponds to one quadratic form; the     */
      /* entries in column 0 and 1 correspond to A and B, respectively.      */
      /* Only forms with B >= 0 are stored, since the others simply yield    */
      /* complex conjugate eta values.                                       */
   int h1, h2, h12, h;
      /* the numbers of ambiguous reduced forms, pairs of non-ambiguous      */
      /* reduced forms, h1+h2 and the class number                           */
} cm_classgroup_t;

typedef struct {
   cm_modular_t m;
   cm_classgroup_t cl;
   mpfr_t root;
      /* sqrt (-cl.d); */
   mpfr_t sqrt2_over2, sqrt2_over4;
   mpc_t *eta_value;
      /* contains the values of eta with respect to the entries of           */
      /* cl.form with the same row index. So only cl.h12 values are stored.  */
} cm_modclass_t;

typedef struct {
   char invariant;
      /* a constant describing which invariant is actually used                 */
   int field;
      /* a constant describing whether we are working over the real or the      */
      /* complex numbers                                                        */
   int p;
      /* some parameter of the class invariant                                  */
   int_cl_t d;
      /* the discriminant                                                       */
   int h, h1, h2, h12;
      /* the class number, the number of forms leading to one or two            */
      /* conjugates, and h1+h2                                                  */
   int minpoly_deg;
      /* the degree of the minimal polynomial; usually h, always h1 + 2 * h2    */
   mpz_t *minpoly;
      /* real part of the minimal polynomial of the function over Q             */
   mpz_t *minpoly_complex;
      /* Only meaningful in the complex case; then                              */
      /* the minimal polynomial is decomposed into two parts over the integral  */
      /* basis [1, sqrt (D)] resp. [1, (1 + sqrt (D))/2]; the first part is in  */
      /* minpoly from "classinvariant", the second one in this variable.        */
} cm_class_t;


#if defined (__cplusplus)
extern "C" {
#endif

/* functions depending on NTL                                                */
extern void cm_ntl_find_factor (mpz_t *res, mpz_t *f, int f_deg, int factor_deg,
   mpz_t p, bool verbose);
extern void cm_ntl_find_root_monic (mpz_t root, mpz_t *f, int deg, mpz_t p,
   bool verbose);
extern mpz_t* cm_ntl_find_roots (mpz_t *f, int deg, mpz_t p, int *no);


/* functions for classgroups of imaginary-quadratic number fields */

extern void cm_classgroup_init (cm_classgroup_t *cl, int_cl_t disc,
   bool checkpoints, bool verbose);
extern void cm_classgroup_clear (cm_classgroup_t *cl);

extern uint_cl_t cm_classgroup_mod (int_cl_t a, uint_cl_t p);
extern int_cl_t cm_classgroup_gcd (int_cl_t a, int_cl_t b);
extern int cm_classgroup_kronecker (int_cl_t a, int_cl_t b);

extern void cm_classgroup_factor (int_cl_t d,
      uint_cl_t *factors, unsigned int *exponents);
extern int_cl_t cm_classgroup_fundamental_discriminant (int_cl_t d);
extern int cm_classgroup_h (int *h1, int *h2, int_cl_t d);

extern void cm_classgroup_reduce (int_cl_t *a, int_cl_t *b, int_cl_t d);
extern void cm_classgroup_compose (int_cl_t *a, int_cl_t *b,
   int_cl_t a1, int_cl_t b1, int_cl_t a2, int_cl_t b2, int_cl_t d);


/* functions for evaluating modular functions at quadratic integers via
   precomputations */

extern void cm_modclass_init (cm_modclass_t *mc, cm_classgroup_t cl,
   mp_prec_t prec, bool checkpoints, bool verbose);
extern void cm_modclass_clear (cm_modclass_t *mc);

extern void cm_modclass_eta_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_f_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_f1_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_j_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_gamma2_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b);
extern void cm_modclass_gamma3_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b);


/* functions for class polynomials */
extern void cm_class_init (cm_class_t *c, int_cl_t disc, char inv,
   bool verbose);
extern void cm_class_clear (cm_class_t *c);
int cm_class_compute_parameter (int_cl_t disc, int inv, bool verbose);

extern void cm_class_write (cm_class_t c);
extern bool cm_class_read (cm_class_t c);

extern void cm_class_compute_minpoly (cm_class_t c, bool checkpoints,
   bool write, bool verbose);
extern mpz_t* cm_class_get_j_mod_P (int_cl_t d, char inv, mpz_t P, int *no,
   const char* modpoldir, bool read, bool verbose);


#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_CLASS_IMPL_H */
