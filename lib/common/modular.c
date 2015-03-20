/*

modular.c - code for evaluating modular functions in floating point arguments

Copyright (C) 2009, 2010, 2015 Andreas Enge

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

#include "cm_common-impl.h"

static void modular_fundamental_matrix (csrcptr z,cm_matrix_t *M);
static void modular_fundamental_domain_matrix (ctype z, cm_matrix_t *M);
static void atkin_eval (cm_modular_t m, ctype rop, ctype op,
   unsigned long int l);

/*****************************************************************************/
/*                                                                           */
/* internal functions                                                        */
/*                                                                           */
/*****************************************************************************/

static void modular_fundamental_matrix (csrcptr z,cm_matrix_t *M)
   /* computes the transformation matrix M such that Mz is in the            */
   /* fundamental domain                                                     */

{
   bool ok = false;
   long int tmp_int;
   ftype tmp_fr;
   ctype local_z;

   finit (tmp_fr, (mp_exp_t) 100);
   cinit (local_z, (mp_exp_t) 100);

   M->a = 1;
   M->b = 0;
   M->c = 0;
   M->d = 1;

   /* determine the matrix from a low precision approximation of z */
   cset (local_z, z);
   while (!ok) {
      ok = true;
      /* obtain -0.5 <= real part <= 0.5 */
      fround (tmp_fr, crealref (local_z));
      tmp_int = -fget_si (tmp_fr);
      cadd_si (local_z, local_z, tmp_int);
      /* multiply M from the left by T^tmp_int */
      M->a += M->c * tmp_int;
      M->b += M->d * tmp_int;

      cnorm (tmp_fr, local_z);
      if (fcmp_d (tmp_fr, 0.999) < 0) {
         /* apply S */
         cneg (local_z, local_z);
         cui_div (local_z, 1ul, local_z);
         tmp_int = M->a;
         M->a = -M->c;
         M->c = tmp_int;
         tmp_int = M->b;
         M->b = -M->d;
         M->d = tmp_int;
         ok = false;
      }
   }

   /* normalise the matrix */
   if (M->c < 0 || (M->c == 0 && M->d < 0)) {
      M->a = -M->a;
      M->b = -M->b;
      M->c = -M->c;
      M->d = -M->d;
   }

   fclear (tmp_fr);
   cclear (local_z);
}

/*****************************************************************************/

static void modular_fundamental_domain_matrix (ctype z, cm_matrix_t *M)
   /* transforms z into the fundamental domain and returns the inverse       */
   /* transformation matrix M                                                */

{
   long int tmp_int;
   ctype tmp_c1;

   cinit (tmp_c1, cget_prec (z));

   modular_fundamental_matrix (z, M);

   /* apply the matrix to z */
   cmul_si (tmp_c1, z, M->a);
   cadd_si (tmp_c1, tmp_c1, M->b);
   cmul_si (z, z, M->c);
   cadd_si (z, z, M->d);
   cdiv (z, tmp_c1, z);

   /* invert and normalize the matrix */
   tmp_int = M->a;
   M->a = M->d;
   M->d = tmp_int;
   M->b = -M->b;
   M->c = -M->c;
   if (M->c < 0 || (M->c == 0 && M->d < 0)) {
      M->a = -M->a;
      M->b = -M->b;
      M->c = -M->c;
      M->d = -M->d;
   }

   cclear (tmp_c1);
}

/*****************************************************************************/

void cm_modular_fundamental_domain (ctype z)
   /* transforms z into the fundamental domain                               */

{
   ctype tmp_c1;
   cm_matrix_t M;

   cinit (tmp_c1, cget_prec (z));

   modular_fundamental_matrix (z, &M);

   /* apply the matrix to z */
   cmul_si (tmp_c1, z, M.a);
   cadd_si (tmp_c1, tmp_c1, M.b);
   cmul_si (z, z, M.c);
   cadd_si (z, z, M.d);
   cdiv (z, tmp_c1, z);

   cclear (tmp_c1);
}

/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_modular_init (cm_modular_t *m, fprec_t prec)

{
   int   i;

   m->prec = prec;

   finit (m->pi, prec);
   cinit (m->twopii, prec);
   cinit (m->log_zeta24, prec);
   cinit (m->zeta48inv, prec);

   fconst_pi (m->pi);
   cset_ui_ui (m->twopii, 0ul, 0ul);
   fmul_2ui (cimagref (m->twopii), m->pi, 1ul);
   cset_ui_ui (m->log_zeta24, 0ul, 0ul);
   fdiv_ui (cimagref (m->log_zeta24), m->pi, 12ul);
   cdiv_ui (m->zeta48inv, m->log_zeta24, 2ul);
   cneg (m->zeta48inv, m->zeta48inv);
   cexp (m->zeta48inv, m->zeta48inv);

   cinit (m->zeta24 [0], prec);
   cset_ui_ui (m->zeta24 [0], 1ul, 0ul);
   cinit (m->zeta24 [1], prec);
   cexp (m->zeta24 [1], m->log_zeta24);
   for (i = 2; i < 24; i++)
   {
      cinit (m->zeta24 [i], prec);
      cmul (m->zeta24 [i], m->zeta24 [i-1], m->zeta24 [1]);
   }

   finit (m->sqrt2, prec);
   fsqrt_ui (m->sqrt2, 2ul);

   cm_qdev_init (&(m->eta), prec);
}

/*****************************************************************************/

void cm_modular_clear (cm_modular_t *m)

{
   int i;

   cclear (m->log_zeta24);
   cclear (m->twopii);
   fclear (m->pi);
   cclear (m->zeta48inv);
   for (i = 0; i < 24; i++)
      cclear (m->zeta24 [i]);
   fclear (m->sqrt2);
   cm_qdev_clear (&(m->eta));
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_modular_eta_transform (cm_modular_t m, ctype rop, ctype z,
     cm_matrix_t M)
   /* computes the value of eta (M z) / eta (z), where M is assumed to be    */
   /* normalised                                                             */

{
   long int c1, zeta_exp, lambda;
   int a = M.a % 48, b = M.b % 24, c = M.c % 24, d = M.d % 24;
      /* If we do not reduce, then there might be an overflow in zeta_exp.  */
      /* M.a is reduced modulo 48 since the formula contains (M.a*M.a-1)/2. */

   /* eta (M z) = epsilon (M) sqrt (Mc z + Md) eta (z). */
   if (M.c == 0) {
      c1 = 1;
      lambda = 1;
   }
   else {
      c1 = M.c;
      lambda = 0;
      while (c1 % 2 == 0) {
         c1 /= 2;
         lambda++;
      }
   }
   zeta_exp = a * b + c * (d * (1 - a*a) - a) + 3 * (c1 % 8) * (a - 1)
         + (3 * lambda * (a*a - 1)) / 2;
   if (cm_nt_kronecker (M.a, c1) == -1)
      zeta_exp += 12;
   zeta_exp %= 24;
   if (zeta_exp < 0)
      zeta_exp += 24;

   /* compute sqrt (M.c * z + M.d) */
   cmul_si (rop, z, M.c);
   cadd_si (rop, rop, M.d);
   csqrt (rop, rop);

   cmul (rop, rop, m.zeta24 [zeta_exp]);
}

/*****************************************************************************/

void cm_modular_eta_series (cm_modular_t m, ctype rop, ctype q_24)
   /* evaluates the power series for eta with q_24 standing for the 24th     */
   /* root of q; uses an automatically optimised addition chain.             */
   /* All computations are carried out with the precision of rop, that must  */
   /* have the same precision for its real and its imaginary part.           */

{
   ctype factor;

   cinit (factor, cget_prec (rop));

   if (fzero_p (cimagref (q_24))) {
      /* avoid cpow_ui, since it calls cpow to determine the sign of 0 */
      fpow_ui (crealref (factor), crealref (q_24), 24ul);
      fset_ui (cimagref (factor), 0ul);
   }
   else
      cpow_ui (factor, q_24, 24ul);
   cm_qdev_eval (factor, m.eta, factor);
   cmul (rop, q_24, factor);

   cclear (factor);
}

/*****************************************************************************/

void cm_modular_eta_series_fr (cm_modular_t m, ftype rop, ftype q_24)
      /* evaluates the power series for eta with q_24 standing for the 24th     */
      /* root of q; uses an automatically optimised addition chain.             */
      /* All computations are carried out with the precision of rop.            */

{
   ftype factor;

   finit (factor, fget_prec (rop));

   fpow_ui (factor, q_24, 24ul);
   cm_qdev_eval_fr (factor, m.eta, factor);
   fmul (rop, q_24, factor);

   fclear (factor);
}

/*****************************************************************************/

void cm_modular_eta_eval (cm_modular_t m, ctype rop, ctype op)
   /* evaluates the eta function at op */

{
   cm_matrix_t M;

   ctype q24, op_local;
   /* Trick: The transformation into the fundamental domain is carried out   */
   /* with the precision of op, the series evaluation with the precision of  */
   /* rop. In this way, the first computation, which has a tendency to lose  */
   /* digits, can be carried out at a higher precision.                      */
   cinit (q24, cget_prec (rop));
   cinit (op_local, cget_prec (op));

   cset (op_local, op);
   modular_fundamental_domain_matrix (op_local, &M);
   cm_modular_eta_transform (m, rop, op_local, M);
   /* workaround to efficiently handle almost real arguments; here, cexp  */
   /* cannot be improved, since the almost zero imaginary part does have an  */
   /* influence on the imaginary part of the result.                         */
   if (!fzero_p (crealref (op_local)) &&
      fget_exp (crealref (op_local))
         < - 0.8 * ((double) fget_prec (crealref (op_local))))
   {
      cm_modular_eta_eval_fr (m, crealref (q24), cimagref (op_local));
      cmul_fr (rop, rop, crealref (q24));
   }
   else
   {
      cmul (q24, m.log_zeta24, op_local);
      cexp (q24, q24);
      cm_modular_eta_series (m, q24, q24);
      cmul (rop, rop, q24);
   }

   cclear (q24);
   cclear (op_local);
}

/*****************************************************************************/

void cm_modular_eta_eval_fr (cm_modular_t m, ftype rop, ftype op)
   /* evaluates the eta function at the purely imaginary argument op*I,      */
   /* of course without transforming into the fundamental domain             */

{
   ftype q24;

   finit (q24, fget_prec (rop));

   fmul (q24, op, cimagref (m.log_zeta24));
   fneg (q24, q24);
   fexp (q24, q24);
   cm_modular_eta_series_fr (m, rop, q24);

   fclear (q24);
}

/*****************************************************************************/

static void atkin_eval (cm_modular_t m, ctype rop, ctype op,
   unsigned long int l)
   /* evaluates eta (z)*eta (lz) */

{
   ctype tmp, rop_local, op_local;

   cinit (tmp, cget_prec (rop));
   cinit (rop_local, cget_prec (rop));
   cinit (op_local, cget_prec (op));

   cmul_ui (op_local, op, l);
   cm_modular_eta_eval (m, rop_local, op_local);
   cm_modular_eta_eval (m, tmp, op);
   cmul (rop, rop_local, tmp);

   cclear (tmp);
   cclear (rop_local);
   cclear (op_local);
}

/*****************************************************************************/

void cm_modular_atkinhecke_eval (cm_modular_t m, ctype rop, ctype op,
   unsigned long int l, unsigned long int r)
   /* evaluates the quotient of the r-th Hecke operator (for r prime),       */
   /* applied to eta (z)*eta (lz), and the function itself, in the argument  */
   /* -1/z (to obtain a function for Gamma^0 (l) instead of Gamma_0 (l))     */
   /* Expresses the numerator as a in transformed arguments.                 */

{
   ctype Mz, tmp, op_local;
   ctype rop_local;
   unsigned long int i;

   cinit (op_local, cget_prec (op));
   cinit (Mz, cget_prec (op));
   cinit (tmp, cget_prec (rop));
   cinit (rop_local, cget_prec (rop));

   cui_div (op_local, 1ul, op);
   cneg (op_local, op_local);
   cset_ui_ui (rop_local, 0ul, 0ul);
   for (i = 0; i < r; i++) {
      cadd_ui (Mz, op_local, 24*i);
      cdiv_ui (Mz, Mz, r);
      atkin_eval (m, tmp, Mz, l);
      cadd (rop_local, rop_local, tmp);
   }
   cdiv_ui (rop_local, rop_local, r);
   cmul_ui (Mz, op_local, r);
   atkin_eval (m, tmp, Mz, l);
   cadd (rop_local, rop_local, tmp);

   atkin_eval (m, tmp, op_local, l);
   cdiv (rop, rop_local, tmp);

   cclear (op_local);
   cclear (Mz);
   cclear (tmp);
   cclear (rop_local);
}

/*****************************************************************************/

void cm_modular_atkinhecke_level_eval (cm_modular_t m, ctype rop, ctype op,
   unsigned long int l)
   /* evaluates Atkin's optimised function for Gamma^0^* (l) */

{
   if (l == 47) {
      cm_modular_atkinhecke_eval (m, rop, op, 47, 17);
      cneg (rop, rop);
   }
   else if (l == 59) {
      ctype z, tmp;
      cinit (z, cget_prec (op));
      cinit (tmp, cget_prec (rop));
      cset (z, op);
      cm_modular_atkinhecke_eval (m, rop, z, 59, 5);
      cm_modular_atkinhecke_eval (m, tmp, z, 59, 29);
      cadd (rop, rop, tmp);
      cadd_ui (rop, rop, 1ul);
      cclear (z);
      cclear (tmp);
   }
   else if (l == 71) {
      ctype z, tmp;
      cinit (z, cget_prec (op));
      cinit (tmp, cget_prec (rop));
      cset (z, op);
      cm_modular_atkinhecke_eval (m, rop, z, 71, 5);
      cm_modular_atkinhecke_eval (m, tmp, z, 71, 29);
      cadd (rop, rop, tmp);
      cadd_ui (rop, rop, 1ul);
      cclear (z);
      cclear (tmp);
   }
   else if (l == 131) {
      cm_modular_atkinhecke_eval (m, rop, op, 131, 61);
      cadd_ui (rop, rop, 1ul);
   }
   else {
      printf ("*** Called cm_modular_atkinhecke_level_eval with level %li, "
         "for which the optimal Atkin invariant is not implemented.\n", l);
      exit (1);
   }

}

/*****************************************************************************/
