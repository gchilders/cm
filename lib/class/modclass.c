/*

modclass.c - code for evaluating modular functions in quadratic arguments

Copyright (C) 2009, 2010, 2011, 2015 Andreas Enge

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

#include "cm_class-impl.h"

#define mpz_sub_si(c, a, b) \
   (b >= 0) ? mpz_sub_ui (c, a, (unsigned long int) b) \
            : mpz_add_ui (c, a, (unsigned long int) (-b))

static void cm_modclass_cset_quadratic (ctype rop,
   int_cl_t a, int_cl_t b, ftype root);
static void write_q24eta (ctype *q24eta, cm_classgroup_t cl,
   fprec_t prec, char *type);
static bool read_q24eta (ctype *q24eta, cm_classgroup_t cl, fprec_t prec,
   char *type);
static void compute_q24 (cm_modular_t m, cm_classgroup_t cl, ftype root,
   ctype *q24, bool verbose);
static void compute_eta (cm_modular_t m, cm_classgroup_t cl, ftype root,
   ctype *eta, bool fem, bool checkpoints, bool verbose);
static void multieta_eval_quad_rec (cm_modclass_t mc, ctype rop_num,
   ctype rop_den, int_cl_t a, int_cl_t b, int *p);

/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

#define CM_FEM false
void cm_modclass_init (cm_modclass_t *mc, cm_classgroup_t cl,
   cm_classgroup_t cl2, fprec_t prec, bool checkpoints, bool verbose)
   /* cl2.d may be 0, in which case mc->eta2 is not computed                 */

{
   int i;
   mpz_t tmp_z;

   mpz_init (tmp_z);

   cm_modular_init (&(mc->m), prec);
   mc->cl = cl;
   mc->cl2 = cl2;

   finit (mc->root, prec);
   cm_classgroup_mpz_set_icl (tmp_z, -cl.d);
   fset_z (mc->root, tmp_z);
   fsqrt (mc->root, mc->root);
   finit (mc->sqrt2_over2, prec);
   fsqrt_ui (mc->sqrt2_over2, 2ul);
   fdiv_2ui (mc->sqrt2_over2, mc->sqrt2_over2, 1ul);
   finit (mc->sqrt2_over4, prec);
   fdiv_2ui (mc->sqrt2_over4, mc->sqrt2_over2, 1ul);

   mc->eta = (ctype *) malloc (mc->cl.h12 * sizeof (ctype));
   for (i = 0; i < mc->cl.h12; i++)
      cinit (mc->eta [i], prec);
   if (!checkpoints || !read_q24eta (mc->eta, mc->cl, prec, "eta"))
      compute_eta (mc->m, mc->cl, mc->root, mc->eta, CM_FEM, checkpoints, verbose);

   if (cl2.d != 0) {
      finit (mc->root2, prec);
      cm_classgroup_mpz_set_icl (tmp_z, -cl2.d);
      fset_z (mc->root2, tmp_z);
      fsqrt (mc->root2, mc->root2);
      mc->eta2 = (ctype *) malloc (mc->cl2.h12 * sizeof (ctype));
      for (i = 0; i < mc->cl2.h12; i++)
         cinit (mc->eta2 [i], prec);
      if (!checkpoints || !read_q24eta (mc->eta2, mc->cl2, prec, "eta"))
         compute_eta (mc->m, mc->cl2, mc->root2, mc->eta2, CM_FEM, checkpoints,
            verbose);
   }

   mpz_clear (tmp_z);
}

/*****************************************************************************/

void cm_modclass_clear (cm_modclass_t *mc)

{
   int i;

   fclear (mc->root);
   fclear (mc->sqrt2_over2);
   fclear (mc->sqrt2_over4);
   for (i = 0; i < mc->cl.h12; i++)
      cclear (mc->eta [i]);
   free (mc->eta);
   if (mc->cl2.d != 0) {
      fclear (mc->root2);
      for (i = 0; i < mc->cl2.h12; i++)
         cclear (mc->eta2 [i]);
      free (mc->eta2);
   }

   cm_modular_clear (&(mc->m));
}

/*****************************************************************************/
/*                                                                           */
/* file handling functions                                                   */
/*                                                                           */
/*****************************************************************************/

static void write_q24eta (ctype *q24eta, cm_classgroup_t cl,
   fprec_t prec, char *type)
   /* writes the values of q24eta to the file                                */
   /* CLASS_TMPDIR + "/tmp_" + -cl.d + "_" + prec + "_" + type + ".dat"      */
   /* type should be one of "q24" or "eta"                                   */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%d_%s.dat", CM_CLASS_TMPDIR, -cl.d,
      (int) prec, type);

   if (!cm_file_open_write (&f, filename))
      exit (1);

   for (i = 0; i < cl.h12; i++) {
      cout_str (f, 16, 0, q24eta [i]);
      fprintf (f, "\n");
   }

   cm_file_close (f);
}

/*****************************************************************************/

static bool read_q24eta (ctype *q24eta, cm_classgroup_t cl, fprec_t prec,
   char *type)
   /* reads the values of q24eta from the file                               */
   /* CLASS_TMPDIR + "/tmp_" + -cl.d + "_" + prec + "_" + type + ".dat"      */
   /* type should be one of "q24" or "eta"                                   */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%d_%s.dat", CM_CLASS_TMPDIR, -cl.d,
            (int) prec, type);

   if (!cm_file_open_read (&f, filename))
      return false;

   for (i = 0; i < cl.h12; i++) {
      cinp_str (q24eta [i], f, NULL, 16);
   }

   cm_file_close (f);

   return true;
}

/*****************************************************************************/
/*                                                                           */
/* functions for precomputations                                             */
/*                                                                           */
/*****************************************************************************/

static void compute_q24 (cm_modular_t m, cm_classgroup_t cl, ftype root,
   ctype *q24, bool verbose)
   /* Computes the q^(1/24) for all forms of cl; root contains sqrt(-d).     */
{
   int i, j, tmp_int;
   int_cl_t tmp_cli;
   ftype Pi24, Pi24_root, tmp;
   ftype *q_real;
   unsigned long int *A_red, *B_red, *order;
   cm_timer  clock2, clock3;
   int    counter1, counter2;

   finit (Pi24, m.prec);
   finit (Pi24_root, m.prec);
   finit (tmp, m.prec);

   fdiv_ui (Pi24, m.pi, 24ul);
   fmul (Pi24_root, root, Pi24);

   A_red = (unsigned long int *) malloc (cl.h12 * sizeof (unsigned long int));
   B_red = (unsigned long int *) malloc (cl.h12 * sizeof (unsigned long int));
   order = (unsigned long int *) malloc (cl.h12 * sizeof (unsigned long int));
   q_real = (ftype *) malloc (cl.h12 * sizeof (ftype));

   for (i = 0; i < cl.h12; i++)
      finit (q_real [i], m.prec);

   cm_timer_start (clock2);
   /* precompute the absolute values of q^{1/24}, i.e. the */
   /* exp (- pi/24 * \sqrt (|d|) / A) for the occurring A  */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   for (i = cl.h12 - 1; i >= 0; i--) {
      /* Check whether the current A is a divisor of a previous one; */
      /* if so, take the previous one which is closest.              */
      /* Notice that the A's are in increasing order.                */
      for (j = i+1;
           j < cl.h12 && cl.form [j].a % cl.form [i].a != 0; j++);
      if (j < cl.h12) {
         if (cl.form [i].a == cl.form [j].a)
            fset (q_real [i], q_real [j]);
         else {
            counter1++;
            fpow_ui (q_real [i], q_real [j],
                         cl.form [j].a / cl.form [i].a);
         }
      }
      else {
         counter1++;
         counter2++;
         fdiv_ui (q_real [i], Pi24_root, cl.form [i].a);
         fneg (q_real [i], q_real [i]);
         fexp (q_real [i], q_real [i]);
      }
      if (verbose &&i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("\n- Number of distinct A:           %d\n", counter1);
      printf ("- Number of exp:                  %d\n", counter2);
      printf ("- Time for exp:                   %.1f\n", cm_timer_get (clock3));
   }

   /* precompute the arguments of q^{1/24}, i.e. the                     */
   /* exp (- i pi/24 * B / A) for the occurring A, B                     */
   /* Write B / A = B_red / A_red with coprime (B_red, A_red).           */
   /* Then order the forms with respect to increasing A_red, and for the */
   /* same A_red with respect to decreasing B_red                        */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   for (i = 0; i < cl.h12; i++) {
      tmp_cli = cm_classgroup_gcd (cl.form [i].a, cl.form [i].b);
      A_red [i] = cl.form [i].a / tmp_cli;
      B_red [i] = cl.form [i].b / tmp_cli;
      order [i] = i;
   }

   /* sort by insertion */
   for (i = 0; i < cl.h12 - 1; i++) {
      tmp_int = i;
      for (j = tmp_int + 1; j < cl.h12; j++)
         if (A_red [order [j]] < A_red [order [tmp_int]]
             || (   A_red [order [j]] == A_red [order [tmp_int]]
             && B_red [order [j]] > B_red [order [tmp_int]]))
            tmp_int = j;
      j = order [i];
      order [i] = order [tmp_int];
      order [tmp_int] = j;
   }

   /* put q24 [] = exp (- pi/24 i / A_red), i.e. the primitive root  */
   /* of unity                                                       */
   for (i = cl.h12 - 1; i >= 0; i--) {
      /* Check whether the current A is a divisor of a previous one; */
      /* if so, take the previous one which is closest.              */
      /* Notice that the A's are in increasing order.                */
      for (j = i+1;
           j < cl.h12 && A_red [order [j]] % A_red [order [i]] != 0;
           j++);
      if (j < cl.h12) {
         if (A_red [order [i]] == A_red [order [j]])
            cset (q24 [order [i]], q24 [order [j]]);
         else {
            counter1++;
            cpow_ui (q24 [order [i]], q24 [order [j]],
                        A_red [order [j]] / A_red [order [i]]);
         }
      }
      else {
         counter1++;
         counter2++;
         fdiv_ui (tmp, Pi24, A_red [order [i]]);
         fsin_cos (q24 [order [i]]->im, q24 [order [i]]->re,
                       tmp);
      }
      if (verbose && i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("\n- Number of distinct A_red:       %d\n", counter1);
      printf ("- Number of sin/cos:              %d\n", counter2);
      printf ("- Time for sin/cos:               %.1f\n", cm_timer_get (clock3));
   }

   /* now raise to the power B_red */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   for (i = cl.h12 - 1; i >= 0; i--)
   {
      if (B_red [order [i]] == 0)
         cset_ui_ui (q24 [order [i]], 1ul, 0ul);
      else if (B_red [order [i]] > 1) {
         /* Check whether the current B is a multiple of a previous one with */
         /* the same A; if so, take the previous one which is closest.       */
         for (j = i+1;
              j < cl.h12 && A_red [order [j]] == A_red [order [i]]
                    && B_red [order [j]] > 0
                    && B_red [order [i]] % B_red [order [j]] != 0;
              j++);
         if (j < cl.h12 && A_red [order [j]] == A_red [order [i]]
             && B_red [order [j]] > 0) {
            if (B_red [order [i]] == B_red [order [j]])
               cset (q24 [order [i]], q24 [order [j]]);
            else {
               counter1++;
               cpow_ui (q24 [order [i]], q24 [order [j]],
                           B_red [order [i]] / B_red [order [j]]);
            }
         }
         else {
            counter2++;
            cpow_ui (q24 [order[i]], q24 [order [i]],
                        B_red [order [i]]);
         }
      }
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("- Number of partial powers:       %d\n", counter1);
      printf ("- Number of total powers:         %d\n", counter2);
      printf ("- Time for powers:                %.1f\n", cm_timer_get (clock3));
   }

   /* Compute the q^(1/24) in q24 */
   for (i = 0; i < cl.h12; i++)
      cmul_fr (q24 [i], q24 [i], q_real [i]);
   cm_timer_stop (clock2);
   if (verbose)
      printf ("- Time for q^(1/24):              %.1f\n", cm_timer_get (clock2));

   for (i = 0; i < cl.h12; i++)
      fclear (q_real [i]);

   free (A_red);
   free (B_red);
   free (order);
   free (q_real);

   fclear (Pi24);
   fclear (Pi24_root);
   fclear (tmp);
}

/*****************************************************************************/

static void compute_eta (cm_modular_t m, cm_classgroup_t cl, ftype root,
   ctype *eta, bool fem, bool checkpoints, bool verbose)
   /* computes the values of the Dedekind eta function for all reduced forms */
   /* of the class group cl and stores them in eta; the precision is taken   */
   /* from eta [0], and root contains sqrt(-d).                              */
   /* If fem is true, then evaluation via the AGM is activated.              */

{
   int i;
   cm_timer clock1;
   fprec_t prec = cget_prec (eta [0]);
   cm_timer_start (clock1);

   if (!fem) {
      ctype *q24;
      cm_timer clock2;

      q24 = (ctype *) malloc (cl.h12 * sizeof (ctype));
      for (i = 0; i < cl.h12; i++) {
         cinit (q24 [i], prec);
      }

      if (!checkpoints || !read_q24eta (q24, cl, prec, "q24")) {
         compute_q24 (m, cl, root, q24, verbose);
         if (checkpoints)
            write_q24eta (q24, cl, prec, "q24");
      }

      cm_timer_start (clock2);
      for (i = 0; i < cl.h12; i++) {
         cm_modular_eta_series (m, eta [i], q24 [i]);
         if (verbose && i % 200 == 0) {
            printf (".");
            fflush (stdout);
         }
      }
      cm_timer_stop (clock2);
      if (verbose)
         printf ("\n- Time for series:                %.1f",
                  cm_timer_get (clock2));
      cm_timer_stop (clock2);
      if (checkpoints)
         write_q24eta (eta, cl, prec, "eta");

      for (i = 0; i < cl.h12; i++)
         cclear (q24 [i]);
      free (q24);
   }
   else {
      ctype tau;

      cinit (tau, prec);

      for (i = 0; i < cl.h12; i++) {
         cm_modclass_cset_quadratic (tau, cl.form [i].a, cl.form [i].b,
            root);
         cm_fem_eta_eval (m, eta [i], tau);
         if (verbose && i % 200 == 0) {
            printf (".");
            fflush (stdout);
         }
      }

      cclear (tau);
   }

   cm_timer_stop (clock1);
   if (verbose)
      printf ("\n- Time for eta:                   %.1f\n",
              cm_timer_get (clock1));
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
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_modclass_eta_eval_quad (ctype rop, cm_modular_t m, cm_classgroup_t cl,
   ctype *eta, int_cl_t a, int_cl_t b, ftype root)
   /* evaluates the eta function at the quadratic integer                    */
   /* (b + sqrt (cl.d)) / (2a)                                               */
   /* root needs to be set to sqrt (|d|).                                    */
   /* eta is a list of precomputed values, indexed in the same way as the    */
   /* list of reduced quadratic forms in cl.                                 */

{
   int_cl_t a_local, b_local;
   int      i, sign;
   cm_matrix_t M;
   ctype    tmp;

   a_local = a;
   b_local = b;
   cm_modclass_fundamental_domain_quad (cl.d, &a_local, &b_local, &M);
   cm_modclass_cset_quadratic (rop, a_local, b_local, root);
   cm_modular_eta_transform (m, rop, rop, M);

   /* look up the eta value */
   i = 0;
   if (b_local < 0) {
      b_local = -b_local;
      sign = -1;
   }
   else
      sign = 1;
   while ((i < cl.h12)
     && (cl.form [i].a != a_local || cl.form [i].b != b_local))
      i++;
   if (i == cl.h12) {
      /* eta value not found, compute it. May happen when the level of the   */
      /* modular function and the conductor have a common factor. This case  */
      /* is rare, and the following computations are not optimised.          */
      printf ("Q");
      cinit (tmp, cget_prec (rop));
      if (sign == 1)
         cm_modclass_cset_quadratic (tmp, a_local, b_local, root);
      else
         cm_modclass_cset_quadratic (tmp, a_local, -b_local, root);
      cm_modular_eta_eval (m, tmp, tmp);
      cmul (rop, rop, tmp);
      cclear (tmp);
   }
   else {
      if (sign == 1)
         cmul (rop, rop, eta [i]);
      else {
         cconj (eta [i], eta [i]);
         cmul (rop, rop, eta [i]);
         cconj (eta [i], eta [i]);
      }
   }
}

/*****************************************************************************/

void cm_modclass_f_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)
   /* evaluates the Weber f-function at the quadratic integer                */
   /* (b + sqrt (d)) / (2*a)                                                 */

{
   ctype tmp;
   int_cl_t c;

   cinit (tmp, cget_prec (rop));

   cm_modclass_eta_eval_quad (tmp, mc.m, mc.cl, mc.eta, a, b, mc.root);

   /* The argument (z+1)/2 of the numerator corresponds to the quadratic     */
   /* form [2*a, b+2*a, (a+b+c)/2]. Here, (a+b+c)/2 need not be integral;    */
   /* if it is, the form need not be primitive any more, but may have a      */
   /* common divisor 2.                                                      */
   c = a + b + cm_classgroup_compute_c (a, b, mc.cl.d);
   if (c % 2 == 0)
      if (b % 2 != 0 || c % 4 != 0)
         cm_modclass_eta_eval_quad (rop, mc.m, mc.cl, mc.eta,
            2*a, b + 2*a, mc.root);
      else {
         assert (mc.cl2.d == mc.cl.d / 4);
         cm_modclass_eta_eval_quad (rop, mc.m, mc.cl2, mc.eta2, a, b/2 + a, mc.root2);
      }
   else {
      ctype z;
      cinit (z, cget_prec (rop));
      cm_modclass_cset_quadratic (z, a, b, mc.root);
      cadd_ui (rop, z, 1ul);
      cdiv_ui (rop, rop, 2ul);
      cm_modular_eta_eval (mc.m, rop, rop);
      cclear (z);
   }

   cdiv (rop, rop, tmp);
   cmul (rop, rop, mc.m.zeta48inv);

   cclear (tmp);
}

/*****************************************************************************/

void cm_modclass_f1_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)
   /* evaluates the Weber f1-function at the quadratic integer               */
   /* (b + sqrt (d)) / (2*a); the same comments as for f apply.              */

{
   ctype tmp;
   int_cl_t c;

   cinit (tmp, cget_prec (rop));

   cm_modclass_eta_eval_quad (tmp, mc.m, mc.cl, mc.eta, a, b, mc.root);
   /* The argument z/2 of the numerator corresponds to the quadratic form    */
   /* [2*a, b, c/2]. Here, c/2 need not be integral; if it is, the form need */
   /* not be primitive any more, but may have a common divisor 2.            */
   c = cm_classgroup_compute_c (a, b, mc.cl.d);
   if (c % 2 == 0)
      if (b % 2 != 0 || c % 4 != 0)
         cm_modclass_eta_eval_quad (rop, mc.m, mc.cl, mc.eta, 2*a, b, mc.root);
      else {
         assert (mc.cl2.d == mc.cl.d / 4);
         cm_modclass_eta_eval_quad (rop, mc.m, mc.cl2, mc.eta2, a, b/2, mc.root2);
      }
   else {
      ctype z;
      cinit (z, cget_prec (rop));
      cm_modclass_cset_quadratic (z, a, b, mc.root);
      cdiv_ui (rop, z, 2ul);
      cm_modular_eta_eval (mc.m, rop, rop);
      cclear (z);
   }
   cdiv (rop, rop, tmp);

   cclear (tmp);
}

/*****************************************************************************/

void cm_modclass_gamma2_eval_quad (cm_modclass_t mc, ctype rop,
   int_cl_t a, int_cl_t b)

{
   ctype f;

   cinit (f, cget_prec (rop));

   cm_modclass_f1_eval_quad (mc, f, a, b);
   cpow_ui (f, f, 8ul);
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

   cm_modclass_f_eval_quad (mc, f, a, b);
   cpow_ui (f, f, 8ul);
   cm_modclass_f1_eval_quad (mc, rop, a, b);
   cpow_ui (rop, rop, 8ul);

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
   ctype rop_den, int_cl_t a, int_cl_t b, int *p)
   /* Evaluates a multiple eta quotient, whose transformation degress are    */
   /* given by the numbers in p, which is an array terminated by 0; the      */
   /* result is given by rop_num / rop_den. This approach replaces complex   */
   /* divisions by faster multiplications.                                   */
   /* Assumes that p contains at least one number.                           */

{
   if (p [1] == 0) {
      /* simple eta quotient */
      cm_modclass_eta_eval_quad (rop_num, mc.m, mc.cl, mc.eta, a * p [0], b,
         mc.root);
      cm_modclass_eta_eval_quad (rop_den, mc.m, mc.cl, mc.eta, a, b, mc.root);
   }
   else if (p [2] == 0 && p [0] == p [1]) {
      /* special, faster code for double eta quotients with twice the same   */
      /* transformation degree                                               */
      cm_modclass_eta_eval_quad (rop_den, mc.m, mc.cl, mc.eta, a, b, mc.root);
      cm_modclass_eta_eval_quad (rop_num, mc.m, mc.cl, mc.eta,
         a * p [0] * p [0], b, mc.root);
      cmul (rop_den, rop_den, rop_num);
      cm_modclass_eta_eval_quad (rop_num, mc.m, mc.cl, mc.eta, a * p [0], b,
         mc.root);
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
   int_cl_t a, int_cl_t b, int *p, int e)
   /* Evaluates a multiple eta quotient, whose transformation degress are    */
   /* given by the numbers in p, which is an array terminated by 0; the      */
   /* quotient is additionally raised to the power e.                        */
   /* Assumes that p contains at least one number.                           */

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
