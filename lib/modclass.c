/*

modclass.c - code for evaluating modular functions in quadratic arguments

Copyright (C) 2009, 2010, 2011, 2015, 2016, 2018, 2021 Andreas Enge

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

#include "cm-impl.h"

#define mpz_sub_si(c, a, b) \
   (b >= 0) ? mpz_sub_ui (c, a, (unsigned long int) b) \
            : mpz_add_ui (c, a, (unsigned long int) (-b))

static int form_cmp (const void *f1, const void *f2);
static void cm_modclass_cset_quadratic (ctype rop,
   int_cl_t a, int_cl_t b, ftype root);
static int cm_modclass_eta_transform_eval_quad (ctype rop, long int *e,
   ctype czplusd, cm_modclass_t mc, int_cl_t a, int_cl_t b);
static void compute_q24 (cm_modclass_t mc, ctype *q24, bool verbose);
static void compute_eta (cm_modclass_t mc, bool verbose);
static void multieta_eval_quad_rec (cm_modclass_t mc, ctype rop_num,
   ctype rop_den, int_cl_t a, int_cl_t b, const int *p);

/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_modclass_init (cm_modclass_t *mc, int_cl_t d, fprec_t prec,
   bool verbose)

{
   int i;
   mpz_t tmp_z;
   cm_classgroup_t cl;

   mpz_init (tmp_z);

   cm_modular_init (&(mc->m), prec);

   mc->d = d;
   finit (mc->root, prec);
   cm_classgroup_mpz_set_icl (tmp_z, -d);
   fset_z (mc->root, tmp_z);
   fsqrt (mc->root, mc->root);

   /* Compute the reduced forms of discriminant d with non-negative b. */
   cm_classgroup_init (&cl, d, false);
   mc->form = (cm_form_t *) malloc (cl.h * sizeof (cm_form_t));
   mc->eta = (ctype *) malloc (cl.h * sizeof (ctype));
   mc->h12 = 0;
   for (i = 0; i < cl.h; i++)
      if (cl.form [i].b >= 0) {
         mc->form [mc->h12] = cl.form [i];
         cinit (mc->eta [mc->h12], prec);
         mc->h12++;
      }
   mc->form = (cm_form_t *) realloc (mc->form,
      mc->h12 * sizeof (cm_form_t));
   mc->eta = (ctype *) realloc (mc->eta, mc->h12 * sizeof (ctype));
   cm_classgroup_clear (&cl);
   qsort (mc->form, mc->h12, sizeof (cm_form_t), form_cmp);

   compute_eta (*mc, verbose);

   mpz_clear (tmp_z);
}

/*****************************************************************************/

void cm_modclass_clear (cm_modclass_t *mc)

{
   int i;

   fclear (mc->root);
   for (i = 0; i < mc->h12; i++)
      cclear (mc->eta [i]);
   free (mc->eta);
   free (mc->form);

   cm_modular_clear (&(mc->m));
}

/*****************************************************************************/
/*                                                                           */
/* Functions for comparing and sorting quadratic forms.                      */
/*                                                                           */
/*****************************************************************************/

static int form_cmp (const void *f1, const void *f2)
   /* All considered forms have a positive a and a non-negative b; they
      are sorted by increasing a and, for the same a, by decreasing b. */
{
   cm_form_t F1, F2;

   F1 = *((cm_form_t *) f1);
   F2 = *((cm_form_t *) f2);

   if (F1.a < F2.a)
      return -1;
   else if (F1.a > F2.a)
      return +1;
   else if (F1.b > F2.b)
      return -1;
   else if (F1.b < F2.b)
      return +1;
   else
      return 0;
}

/*****************************************************************************/
/*                                                                           */
/* functions for precomputations                                             */
/*                                                                           */
/*****************************************************************************/

static void compute_q24 (cm_modclass_t mc, ctype *q24, bool verbose)
   /* Compute the q^(1/24) for all forms in mc and return them in q24. */
{
   int i, j;
   ftype Pi24, Pi24_root, tmp;
   ftype *q_real;
   cm_timer clock2, clock3;
   int counter1, counter2;

   finit (Pi24, mc.m.prec);
   finit (Pi24_root, mc.m.prec);
   finit (tmp, mc.m.prec);

   fdiv_ui (Pi24, mc.m.pi, 24ul);
   fmul (Pi24_root, mc.root, Pi24);
   fneg (Pi24_root, Pi24_root);

   cm_timer_start (clock2);
   /* Compute in q_real the absolute values of q^{1/24}, i.e. the
      exp (- pi/24 * \sqrt (|d|) / a) for the occurring a,
      and at the same time in q24 the roots of unity
      exp (pi/24 * I / a). */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   q_real = (ftype *) malloc (mc.h12 * sizeof (ftype));
   for (i = 0; i < mc.h12; i++)
      finit (q_real [i], mc.m.prec);
   for (i = mc.h12 - 1; i >= 0; i--) {
      /* Check whether the current a is a divisor of a previous one.
         Then the q-values can be computed by an exponentiation. The
         exponent is minimised by choosing the closest a that is divided. */
      for (j = i+1;
           j < mc.h12 && mc.form [j].a % mc.form [i].a != 0;
           j++);
      if (j < mc.h12) {
         if (mc.form [i].a == mc.form [j].a) {
            fset (q_real [i], q_real [j]);
            cset (q24 [i], q24 [j]);
         }
         else {
            counter1++;
            fpow_ui (q_real [i], q_real [j], mc.form [j].a / mc.form [i].a);
            cpow_ui (q24 [i], q24 [j], mc.form [j].a / mc.form [i].a);
         }
      }
      else {
         counter1++;
         counter2++;
         fdiv_ui (q_real [i], Pi24_root, mc.form [i].a);
         fexp (q_real [i], q_real [i]);
         fdiv_ui (tmp, Pi24, mc.form [i].a);
         fsin_cos (q24 [i]->im, q24 [i]->re, tmp);
      }
   }
   if (verbose && i % 200 == 0) {
      printf (".");
      fflush (stdout);
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("\n- Number of distinct A:  %d\n", counter1);
      printf ("- Number of exp/sin_cos: %d\n", counter2);
      printf ("- Time for exp/sin_cos:  %.1f\n",
         cm_timer_get (clock3));
   }

   /* Raise the roots of unity in q24 to the powers b. */
   cm_timer_start (clock3);
   for (i = 0; i < mc.h12; i++)
      if (mc.form [i].b == 0)
         cset_ui_ui (q24 [i], 1ul, 0ul);
      else if (mc.form [i].b > 1)
         cpow_ui (q24 [i], q24 [i], (unsigned long int) mc.form [i].b);
   cm_timer_stop (clock3);
   printf ("- Time for B powers:     %.1f\n", cm_timer_get (clock3));

   /* Compute the q^(1/24) in q24 */
   for (i = 0; i < mc.h12; i++)
      cmul_fr (q24 [i], q24 [i], q_real [i]);

   cm_timer_stop (clock2);
   if (verbose)
      printf ("- Time for q^(1/24):     %.1f\n", cm_timer_get (clock2));

   for (i = 0; i < mc.h12; i++)
      fclear (q_real [i]);
   free (q_real);

   fclear (Pi24);
   fclear (Pi24_root);
   fclear (tmp);
}

/*****************************************************************************/

static void compute_eta (cm_modclass_t mc, bool verbose)
   /* Compute the values of the Dedekind eta function for all reduced forms
      in mc and store them in mc. */
{
   int i;
   cm_timer clock1;
   fprec_t prec = mc.m.prec;

   ctype *q24;
   cm_timer clock2;

   cm_timer_start (clock1);
   q24 = (ctype *) malloc (mc.h12 * sizeof (ctype));
   for (i = 0; i < mc.h12; i++)
      cinit (q24 [i], prec);

   compute_q24 (mc, q24, verbose);

   cm_timer_start (clock2);
   for (i = 0; i < mc.h12; i++) {
      cm_modular_eta_series (mc.m, mc.eta [i], q24 [i]);
      if (verbose && i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   cm_timer_stop (clock2);
   if (verbose)
      printf ("\n- Time for series:       %.1f", cm_timer_get (clock2));
   cm_timer_stop (clock2);

   for (i = 0; i < mc.h12; i++)
      cclear (q24 [i]);
   free (q24);

   cm_timer_stop (clock1);
   if (verbose)
      printf ("\n- Time for eta:          %.1f\n", cm_timer_get (clock1));
}

/*****************************************************************************/
/*                                                                           */
/* other internal functions                                                  */
/*                                                                           */
/*****************************************************************************/

static void cm_modclass_cset_quadratic (ctype rop,
   int_cl_t a, int_cl_t b, ftype root)
   /* sets rop to (b + i*root) / (2a)                                      */

{
   fset_si (rop->re, b);
   fset (rop->im, root);
   cdiv_ui (rop, rop, 2*a);
}

/*****************************************************************************/

static void cm_modclass_fundamental_domain_quad (int_cl_t d, int_cl_t *a,
   int_cl_t *b, cm_matrix_t *M)
      /* transforms (b + sqrt (d)) / (2a) into the fundamental domain and   */
      /* returns the inverse transformation matrix M                        */
      /* Unfortunately we have to compute internally with multiprecision;   */
      /* although the final result fits into a long, intermediate results   */
      /* can be quite large.                                                */
{
   bool reduced = false;
   static mpz_t a_local, b_local, b_minus_a, two_a, offset, c;
   static bool first_time = true;
   int_cl_t     tmp_int;

   if (first_time)
   {
      mpz_init (a_local);
      mpz_init (b_local);
      mpz_init (b_minus_a);
      mpz_init (two_a);
      mpz_init (offset);
      mpz_init (c);
      first_time = false;
   }

   M->a = 1;
   M->b = 0;
   M->c = 0;
   M->d = 1;
   cm_classgroup_mpz_set_icl (a_local, *a);
   cm_classgroup_mpz_set_icl (b_local, *b);

   while (!reduced) {
      /* obtain -a < b <= a                                                   */
      /* basically, compute offset as b / (2a) rounded towards the closest    */
      /* integer. If the quotient is exactly between two integers, round down */
      mpz_sub (b_minus_a, b_local, a_local);
      mpz_mul_2exp (two_a, a_local, 1);
      mpz_cdiv_q (offset, b_minus_a, two_a);

      tmp_int = mpz_get_si (offset);
      /* multiply M from the right by T^{tmp_int} */
      M->b += M->a * tmp_int;
      M->d += M->c * tmp_int;

      mpz_mul (offset, offset, two_a);
      mpz_sub (b_local, b_local, offset);

      /* compute c */
      mpz_mul (c, b_local, b_local);
      mpz_sub_si (c, c, d);
      mpz_fdiv_q_2exp (c, c, 2);
      mpz_divexact (c, c, a_local);
      /* if not reduced, invert */
      if (mpz_cmp (a_local, c) < 0 ||
          (mpz_cmp (a_local, c) == 0 && mpz_sgn (b_local) >= 0))
         reduced = true;
      else {
         mpz_set (a_local, c);
         mpz_neg (b_local, b_local);
         tmp_int = M->a;
         M->a = M->b;
         M->b = -tmp_int;
         tmp_int = M->c;
         M->c = M->d;
         M->d = -tmp_int;
         reduced = false;
      }
   }

   *a = mpz_get_si (a_local);
   *b = mpz_get_si (b_local);
   /* normalise the matrix */
   if (M->c < 0 || (M->c == 0 && M->d < 0)) {
      M->a = -M->a;
      M->b = -M->b;
      M->c = -M->c;
      M->d = -M->d;
   }
}

/*****************************************************************************/

static int cm_modclass_eta_transform_eval_quad (ctype rop, long int *e,
   ctype czplusd, cm_modclass_t mc, int_cl_t a, int_cl_t b)
   /* Assume that mc is initialised with the discriminant d and that the
      corresponding eta values have been precomputed.
      Transform the quadratic integer op = (b + sqrt (mc.d)) / (2a) into the
      element op1 in the fundamental domain.
      Return eta (op1) in rop and e and czplusd such that
      eta (op) = zeta_24^e * sqrt (czplusd) * rop.
      This can be used to decide at a higher level what to do with the
      transformation data (for instance, if an even exponent is to be
      applied, the square root can be saved).
      Return 0 if op itself is already in the fundamental domain,
      1 otherwise. */
{
   int_cl_t    a_local, b_local;
   int         transformed, i, sign;
   cm_matrix_t M;

   a_local = a;
   b_local = b;
   cm_modclass_fundamental_domain_quad (mc.d, &a_local, &b_local, &M);
   cm_modclass_cset_quadratic (czplusd, a_local, b_local, mc.root);
   transformed = cm_modular_eta_transform (e, czplusd, czplusd, M);

   /* look up the eta value */
   if (b_local < 0) {
      b_local = -b_local;
      sign = -1;
   }
   else
      sign = 1;
   for (i = 0;
      i < mc.h12 && (mc.form [i].a != a_local || mc.form [i].b != b_local);
      i++);
   if (i == mc.h12) {
      /* eta value not found, compute it. May happen when the level of the   */
      /* modular function and the conductor have a common factor. This case  */
      /* is rare, and the following computations are not optimised.          */
      printf ("Q");
      if (sign == 1)
         cm_modclass_cset_quadratic (rop, a_local, b_local, mc.root);
      else
         cm_modclass_cset_quadratic (rop, a_local, -b_local, mc.root);
      cm_modular_eta_eval (mc.m, rop, rop);
   }
   else {
      if (sign == 1)
         cset (rop, mc.eta [i]);
      else
         cconj (rop, mc.eta [i]);
   }

   return (transformed);
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_modclass_eta_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)
   /* Evaluate the eta function at the quadratic integer
      (b + sqrt (mc.d)) / (2a). */
{
   ctype    czplusd;
   long int e;
   int      transformed;

   cinit (czplusd, cget_prec (rop));
   transformed = cm_modclass_eta_transform_eval_quad (rop, &e, czplusd, mc,
      a, b);
   if (transformed) {
      csqrt (czplusd, czplusd);
      cmul (rop, rop, czplusd);
      cmul (rop, rop, mc.m.zeta24 [e]);
   }
}

/*****************************************************************************/

void cm_modclass_f_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, int e)
   /* Evaluate the Weber f-function, raised to the power e, at the quadratic
     integer (b + sqrt (mc.d)) / (2*a).
     e may be 1 or 2. */

{
   ctype    rop1, czplusd1, czplusd2;
   long int e1, e2;
   int_cl_t c;

   cinit (rop1, cget_prec (rop));
   cinit (czplusd1, cget_prec (rop));
   cinit (czplusd2, cget_prec (rop));

   cm_modclass_eta_transform_eval_quad (rop1, &e1, czplusd1, mc, a, b);

   /* The argument (z+1)/2 of the numerator corresponds to the quadratic
      form [2*a, b+2*a, (a+b+c)/2]. Here, (a+b+c)/2 need not be integral;
      if it is, the form need not be primitive any more, but may have a
      common divisor 2. */
   c = a + b + cm_classgroup_compute_c (a, b, mc.d);
   if (c % 2 == 0 && (b % 2 != 0 || c % 4 != 0))
         cm_modclass_eta_transform_eval_quad (rop, &e2, czplusd2, mc,
            2*a, b + 2*a);
   else {
      ctype z;
      cinit (z, cget_prec (rop));
      cm_modclass_cset_quadratic (z, a, b, mc.root);
      cadd_ui (rop, z, 1ul);
      cdiv_ui (rop, rop, 2ul);
      cm_modular_eta_eval (mc.m, rop, rop);
      cclear (z);
      e2 = 0;
      cset_ui (czplusd2, 1ul);
   }

   if (e == 2) {
      csqr (rop1, rop1);
      csqr (rop, rop);
      e2 = (2 * (e2 + (24 - e1)) + 23) % 24;
         /* Force a positive result; 23 stands for the -1 coming from the
            square of zeta48inv. */
   }
   else /* e == 1 */ {
      cmul (rop, rop, mc.m.zeta48inv);
      csqrt (czplusd1, czplusd1);
      csqrt (czplusd2, czplusd2);
      e2 = (e2 + (24 - e1)) % 24;
   }
   cmul (rop1, rop1, czplusd1);
   cmul (rop, rop, czplusd2);
   if (e2 != 0)
      cmul (rop, rop, mc.m.zeta24 [e2]);

   cdiv (rop, rop, rop1);

   cclear (rop1);
   cclear (czplusd1);
   cclear (czplusd2);
}

/*****************************************************************************/

void cm_modclass_f1_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, int e)
   /* Evaluate the Weber f-function, raised to the power e, at the quadratic
     integer (b + sqrt (mc.d)) / (2*a).
     e may be 1 or 2. */

{
   ctype    rop1, czplusd1, czplusd2;
   long int e1, e2;
   int_cl_t c;

   cinit (rop1, cget_prec (rop));
   cinit (czplusd1, cget_prec (rop));
   cinit (czplusd2, cget_prec (rop));

   cm_modclass_eta_transform_eval_quad (rop1, &e1, czplusd1, mc, a, b);

   /* The argument z/2 of the numerator corresponds to the quadratic form
      [2*a, b, c/2]. Here, c/2 need not be integral; if it is, the form need
      not be primitive any more, but may have a common divisor 2. */
   c = cm_classgroup_compute_c (a, b, mc.d);
   if (c % 2 == 0 && (b % 2 != 0 || c % 4 != 0))
         cm_modclass_eta_transform_eval_quad (rop, &e2, czplusd2, mc,
            2*a, b);
   else {
      ctype z;
      cinit (z, cget_prec (rop));
      cm_modclass_cset_quadratic (z, a, b, mc.root);
      cdiv_ui (rop, z, 2ul);
      cm_modular_eta_eval (mc.m, rop, rop);
      cclear (z);
      e2 = 0;
      cset_ui (czplusd2, 1ul);
   }

   if (e == 2) {
      csqr (rop1, rop1);
      csqr (rop, rop);
      e2 = (2 * (e2 + (24 - e1))) % 24;
   }
   else /* e == 1 */ {
      csqrt (czplusd1, czplusd1);
      csqrt (czplusd2, czplusd2);
      e2 = (e2 + (24 - e1)) % 24;
   }
   cmul (rop1, rop1, czplusd1);
   cmul (rop, rop, czplusd2);
   if (e2 != 0)
      cmul (rop, rop, mc.m.zeta24 [e2]);

   cdiv (rop, rop, rop1);

   cclear (rop1);
   cclear (czplusd1);
   cclear (czplusd2);
}

/*****************************************************************************/

void cm_modclass_gamma2_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)

{
   ctype f;

   cinit (f, cget_prec (rop));

   cm_modclass_f1_eval_quad (mc, f, a, b, 2);
   cpow_ui (f, f, 4ul);
   csqr (rop, f);
   cui_div (f, 16ul, f);
   cadd (rop, rop, f);

   cclear (f);
}

/*****************************************************************************/

void cm_modclass_gamma3_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)
   /* evaluates sqrt (d) * gamma3 */

{
   ctype f, tmp;
   ftype tmp_fr;

   cinit (f, cget_prec (rop));
   cinit (tmp, cget_prec (rop));

   cm_modclass_f_eval_quad (mc, f, a, b, 2);
   cpow_ui (f, f, 4ul);
   cm_modclass_f1_eval_quad (mc, rop, a, b, 2);
   cpow_ui (rop, rop, 4ul);

   cmul_ui (rop, rop, 2ul);
   csub (rop, rop, f);
   cpow_ui (tmp, f, 3ul);
   cadd_ui (tmp, tmp, 8ul);
   cmul (rop, rop, tmp);
   cdiv (rop, rop, f);
   cmul_fr (rop, rop, mc.root);
   /* multiply by i */
   tmp_fr [0] = rop->im [0];
   rop->im [0] = rop->re [0];
   rop->re [0] = tmp_fr [0];
   fneg (rop->re, rop->re);

   cclear (f);
   cclear (tmp);
}

/*****************************************************************************/

void cm_modclass_j_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)

{
   cm_modclass_gamma2_eval_quad (mc, rop, a, b);
   cpow_ui (rop, rop, 3ul);
}

/*****************************************************************************/

static void multieta_eval_quad_rec (cm_modclass_t mc, ctype rop_num,
   ctype rop_den, int_cl_t a, int_cl_t b, const int *p)
   /* Evaluate a multiple eta quotient, whose transformation degrees are
      given by the numbers in p, a 0-terminated array with at least one
      entry.
      The result is given by rop_num / rop_den. This approach replaces
      complex divisions by faster multiplications. */
{
   if (p [1] == 0) {
      /* Simple eta quotient. */
      cm_modclass_eta_eval_quad (mc, rop_num, a * p [0], b);
      cm_modclass_eta_eval_quad (mc, rop_den, a, b);
   }
   else if (p [2] == 0 && p [0] == p [1]) {
      /* Special, faster code for double eta quotients with twice the same
         transformation degree. */
      cm_modclass_eta_eval_quad (mc, rop_den, a, b);
      cm_modclass_eta_eval_quad (mc, rop_num, a * p [0] * p [0], b);
      cmul (rop_den, rop_den, rop_num);
      cm_modclass_eta_eval_quad (mc, rop_num, a * p [0], b);
      csqr (rop_num, rop_num);
   }
   else {
      ctype tmp1, tmp2;
      fprec_t prec = cget_prec (rop_num);

      cinit (tmp1, prec);
      cinit (tmp2, prec);

      multieta_eval_quad_rec (mc, rop_num, tmp1, a, b, p+1);
      multieta_eval_quad_rec (mc, rop_den, tmp2, a * p [0], b, p+1);
      cmul (rop_num, rop_num, tmp2);
      cmul (rop_den, rop_den, tmp1);

      cclear (tmp1);
      cclear (tmp2);
   }
}

/*****************************************************************************/

void cm_modclass_multieta_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, const int *p, int e)
   /* Evaluate a multiple eta quotient, whose transformation degrees are
      given by the numbers in p, a 0-terminated array with at least one
      entry. The quotient is additionally raised to the power e. */
{
   ctype tmp;

   cinit (tmp, cget_prec (rop));

   multieta_eval_quad_rec (mc, rop, tmp, a, b, p);
   cdiv (rop, rop, tmp);
   if (e != 1)
      cpow_ui (rop, rop, e);

   cclear (tmp);
}

/*****************************************************************************/

void cm_modclass_atkinhecke_level_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b, unsigned long int l)

{
   ctype z;

   cinit (z, cget_prec (rop));

   cm_modclass_cset_quadratic (z, a, b, mc.root);
   cm_modular_atkinhecke_level_eval (mc.m, rop, z, l);

   cclear (z);
}

/*****************************************************************************/
