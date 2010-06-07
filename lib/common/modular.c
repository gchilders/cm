/*

modular.c - code for evaluating modular functions in floating point arguments

Copyright (C) 2009, 2010 Andreas Enge

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

static void modular_fundamental_matrix (mpc_srcptr z,cm_matrix_t *M);
static void modular_fundamental_domain_matrix (mpc_t z, cm_matrix_t *M);
static void atkin_eval (cm_modular_t m, mpc_t rop, mpc_t op,
   unsigned long int l);

/*****************************************************************************/
/*                                                                           */
/* internal functions                                                        */
/*                                                                           */
/*****************************************************************************/

static void modular_fundamental_matrix (mpc_srcptr z,cm_matrix_t *M)
   /* computes the transformation matrix M such that Mz is in the            */
   /* fundamental domain                                                     */

{
   bool ok = false;
   long int tmp_int;
   mpfr_t tmp_fr;
   mpc_t local_z;

   mpfr_init2 (tmp_fr, (mp_exp_t) 100);
   mpc_init2 (local_z, (mp_exp_t) 100);

   M->a = 1;
   M->b = 0;
   M->c = 0;
   M->d = 1;

   /* determine the matrix from a low precision approximation of z */
   mpc_set (local_z, z, MPC_RNDNN);
   while (!ok) {
      ok = true;
      /* obtain -0.5 <= real part <= 0.5 */
      mpfr_round (tmp_fr, local_z->re);
      tmp_int = -mpfr_get_si (tmp_fr, GMP_RNDN);
      mpc_add_si (local_z, local_z, tmp_int, MPC_RNDNN);
      /* multiply M from the left by T^tmp_int */
      M->a += M->c * tmp_int;
      M->b += M->d * tmp_int;

      mpc_norm (tmp_fr, local_z, MPC_RNDNN);
      if (mpfr_cmp_d (tmp_fr, 0.999) < 0) {
         /* apply S */
         mpc_neg (local_z, local_z, MPC_RNDNN);
         mpc_ui_div (local_z, 1ul, local_z, MPC_RNDNN);
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

   mpfr_clear (tmp_fr);
   mpc_clear (local_z);
}

/*****************************************************************************/

static void modular_fundamental_domain_matrix (mpc_t z, cm_matrix_t *M)
   /* transforms z into the fundamental domain and returns the inverse       */
   /* transformation matrix M                                                */

{
   long int tmp_int;
   mpc_t tmp_c1;

   mpc_init2 (tmp_c1, mpfr_get_prec (z->re));

   modular_fundamental_matrix (z, M);

   /* apply the matrix to z */
   mpc_mul_si (tmp_c1, z, M->a, MPC_RNDNN);
   mpc_add_si (tmp_c1, tmp_c1, M->b, MPC_RNDNN);
   mpc_mul_si (z, z, M->c, MPC_RNDNN);
   mpc_add_si (z, z, M->d, MPC_RNDNN);
   mpc_div (z, tmp_c1, z, MPC_RNDNN);

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

   mpc_clear (tmp_c1);
}

/*****************************************************************************/

void cm_modular_fundamental_domain (mpc_t z)
   /* transforms z into the fundamental domain                               */

{
   mpc_t tmp_c1;
   cm_matrix_t M;

   mpc_init2 (tmp_c1, mpfr_get_prec (z->re));

   modular_fundamental_matrix (z, &M);

   /* apply the matrix to z */
   mpc_mul_si (tmp_c1, z, M.a, MPC_RNDNN);
   mpc_add_si (tmp_c1, tmp_c1, M.b, MPC_RNDNN);
   mpc_mul_si (z, z, M.c, MPC_RNDNN);
   mpc_add_si (z, z, M.d, MPC_RNDNN);
   mpc_div (z, tmp_c1, z, MPC_RNDNN);

   mpc_clear (tmp_c1);
}

/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_modular_init (cm_modular_t *m, mp_prec_t prec)

{
   int   i;

   m->prec = prec;

   mpfr_init2 (m->pi, prec);
   mpc_init2 (m->twopii, prec);
   mpc_init2 (m->log_zeta24, prec);
   mpc_init2 (m->zeta48inv, prec);

   mpfr_const_pi (m->pi, GMP_RNDN);
   mpc_set_ui_ui (m->twopii, 0ul, 0ul, MPC_RNDNN);
   mpfr_mul_2exp (m->twopii->im, m->pi, 1ul, GMP_RNDN);
   mpc_set_ui_ui (m->log_zeta24, 0ul, 0ul, MPC_RNDNN);
   mpfr_div_ui (m->log_zeta24->im, m->pi, 12ul, GMP_RNDN);
   mpc_div_ui (m->zeta48inv, m->log_zeta24, 2ul, MPC_RNDNN);
   mpc_neg (m->zeta48inv, m->zeta48inv, MPC_RNDNN);
   mpc_exp (m->zeta48inv, m->zeta48inv, MPC_RNDNN);

   mpc_init2 (m->zeta24 [0], prec);
   mpc_set_ui_ui (m->zeta24 [0], 1ul, 0ul, MPC_RNDNN);
   mpc_init2 (m->zeta24 [1], prec);
   mpc_exp (m->zeta24 [1], m->log_zeta24, MPC_RNDNN);
   for (i = 2; i < 24; i++)
   {
      mpc_init2 (m->zeta24 [i], prec);
      mpc_mul (m->zeta24 [i], m->zeta24 [i-1], m->zeta24 [1], MPC_RNDNN);
   }

   mpfr_init2 (m->sqrt2, prec);
   mpfr_sqrt_ui (m->sqrt2, 2ul, GMP_RNDN);

   cm_qdev_init (&(m->eta));
}

/*****************************************************************************/

void cm_modular_clear (cm_modular_t *m)

{
   int i;

   mpc_clear (m->log_zeta24);
   mpc_clear (m->twopii);
   mpfr_clear (m->pi);
   mpc_clear (m->zeta48inv);
   for (i = 0; i < 24; i++)
      mpc_clear (m->zeta24 [i]);
   mpfr_clear (m->sqrt2);
   cm_qdev_clear (&(m->eta));
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_modular_eta_transform (cm_modular_t m, mpc_t rop, mpc_t z,
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
   mpc_mul_si (rop, z, M.c, MPC_RNDNN);
   mpc_add_si (rop, rop, M.d, MPC_RNDNN);
   mpc_sqrt (rop, rop, MPC_RNDNN);

   mpc_mul (rop, rop, m.zeta24 [zeta_exp], MPC_RNDNN);
}

/*****************************************************************************/

void cm_modular_eta_series (cm_modular_t m, mpc_t rop, mpc_t q_24)
   /* evaluates the power series for eta with q_24 standing for the 24th     */
   /* root of q; uses an automatically optimised addition chain.             */
   /* All computations are carried out with the precision of rop, that must  */
   /* have the same precision for its real and its imaginary part.           */

{
   mpc_t factor;

   mpc_init2 (factor, mpfr_get_prec (rop->re));

   if (mpfr_zero_p (mpc_imagref (q_24))) {
      /* avoid mpc_pow_ui, since it calls mpc_pow to determine the sign of 0 */
      mpfr_pow_ui (mpc_realref (factor), mpc_realref (q_24), 24ul, GMP_RNDN);
      mpfr_set_ui (mpc_imagref (factor), 0ul, GMP_RNDN);
   }
   else
      mpc_pow_ui (factor, q_24, 24ul, MPC_RNDNN);
   cm_qdev_eval (factor, m.eta, factor);
   mpc_mul (rop, q_24, factor, MPC_RNDNN);

   mpc_clear (factor);
}

/*****************************************************************************/

void cm_modular_eta_series_fr (cm_modular_t m, mpfr_t rop, mpfr_t q_24)
      /* evaluates the power series for eta with q_24 standing for the 24th     */
      /* root of q; uses an automatically optimised addition chain.             */
      /* All computations are carried out with the precision of rop.            */

{
   mpfr_t factor;

   mpfr_init2 (factor, mpfr_get_prec (rop));

   mpfr_pow_ui (factor, q_24, 24ul, GMP_RNDN);
   cm_qdev_eval_fr (factor, m.eta, factor);
   mpfr_mul (rop, q_24, factor, GMP_RNDN);

   mpfr_clear (factor);
}

/*****************************************************************************/

void cm_modular_eta_eval (cm_modular_t m, mpc_t rop, mpc_t op)
   /* evaluates the eta function at op */

{
   cm_matrix_t M;

   mpc_t q24, op_local;
   /* Trick: The transformation into the fundamental domain is carried out   */
   /* with the precision of op, the series evaluation with the precision of  */
   /* rop. In this way, the first computation, which has a tendency to lose  */
   /* digits, can be carried out at a higher precision.                      */
   mpc_init2 (q24, mpfr_get_prec (rop->re));
   mpc_init2 (op_local, mpfr_get_prec (op->re));

   mpc_set (op_local, op, MPC_RNDNN);
   modular_fundamental_domain_matrix (op_local, &M);
   cm_modular_eta_transform (m, rop, op_local, M);
   /* workaround to efficiently handle almost real arguments; here, mpc_exp  */
   /* cannot be improved, since the almost zero imaginary part does have an  */
   /* influence on the imaginary part of the result.                         */
   if (!mpfr_zero_p (op_local->re) &&
      mpfr_get_exp (op_local->re)
         < - 0.8 * ((double) mpfr_get_prec (op_local->re)))
   {
      cm_modular_eta_eval_fr (m, q24->re, op_local->im);
      mpc_mul_fr (rop, rop, q24->re, MPC_RNDNN);
   }
   else
   {
      mpc_mul (q24, m.log_zeta24, op_local, MPC_RNDNN);
      mpc_exp (q24, q24, MPC_RNDNN);
      cm_modular_eta_series (m, q24, q24);
      mpc_mul (rop, rop, q24, MPC_RNDNN);
   }

   mpc_clear (q24);
   mpc_clear (op_local);
}

/*****************************************************************************/

void cm_modular_eta_eval_fr (cm_modular_t m, mpfr_t rop, mpfr_t op)
   /* evaluates the eta function at the purely imaginary argument op*I,      */
   /* of course without transforming into the fundamental domain             */

{
   mpfr_t q24;

   mpfr_init2 (q24, mpfr_get_prec (rop));

   mpfr_mul (q24, op, m.log_zeta24->im, GMP_RNDN);
   mpfr_neg (q24, q24, GMP_RNDN);
   mpfr_exp (q24, q24, GMP_RNDN);
   cm_modular_eta_series_fr (m, rop, q24);

   mpfr_clear (q24);
}

/*****************************************************************************/

static void atkin_eval (cm_modular_t m, mpc_t rop, mpc_t op,
   unsigned long int l)
   /* evaluates eta (z)*eta (lz) */

{
   mpc_t tmp, rop_local, op_local;

   mpc_init2 (tmp, mpc_get_prec (rop));
   mpc_init2 (rop_local, mpc_get_prec (rop));
   mpc_init2 (op_local, mpc_get_prec (op));

   mpc_mul_ui (op_local, op, l, MPC_RNDNN);
   cm_modular_eta_eval (m, rop_local, op_local);
   cm_modular_eta_eval (m, tmp, op);
   mpc_mul (rop, rop_local, tmp, MPC_RNDNN);

   mpc_clear (tmp);
   mpc_clear (rop_local);
   mpc_clear (op_local);
}

/*****************************************************************************/

void cm_modular_atkinhecke_eval (cm_modular_t m, mpc_t rop, mpc_t op,
   unsigned long int l, unsigned long int r)
   /* evaluates the quotient of the r-th Hecke operator (for r prime),       */
   /* applied to eta (z)*eta (lz), and the function itself, in the argument  */
   /* -1/z (to obtain a function for Gamma^0 (l) instead of Gamma_0 (l))     */
   /* Expresses the numerator as a in transformed arguments.                 */

{
   mpc_t Mz, tmp, op_local;
   mpc_t rop_local;
   unsigned long int i;

   mpc_init2 (op_local, mpc_get_prec (op));
   mpc_init2 (Mz, mpc_get_prec (op));
   mpc_init2 (tmp, mpc_get_prec (rop));
   mpc_init2 (rop_local, mpc_get_prec (rop));

   mpc_ui_div (op_local, 1ul, op, MPC_RNDNN);
   mpc_neg (op_local, op_local, MPC_RNDNN);
   mpc_set_ui_ui (rop_local, 0ul, 0ul, MPC_RNDNN);
   for (i = 0; i < r; i++) {
      mpc_add_ui (Mz, op_local, 24*i, MPC_RNDNN);
      mpc_div_ui (Mz, Mz, r, GMP_RNDN);
      atkin_eval (m, tmp, Mz, l);
      mpc_add (rop_local, rop_local, tmp, MPC_RNDNN);
   }
   mpc_div_ui (rop_local, rop_local, r, MPC_RNDNN);
   mpc_mul_ui (Mz, op_local, r, GMP_RNDN);
   atkin_eval (m, tmp, Mz, l);
   mpc_add (rop_local, rop_local, tmp, MPC_RNDNN);

   atkin_eval (m, tmp, op_local, l);
   mpc_div (rop, rop_local, tmp, MPC_RNDNN);

   mpc_clear (op_local);
   mpc_clear (Mz);
   mpc_clear (tmp);
   mpc_clear (rop_local);
}

/*****************************************************************************/

void cm_modular_atkinhecke_level_eval (cm_modular_t m, mpc_t rop, mpc_t op,
   unsigned long int l)
   /* evaluates Atkin's optimised function for Gamma^0^* (l) */

{
   if (l == 47) {
      cm_modular_atkinhecke_eval (m, rop, op, 47, 17);
      mpc_neg (rop, rop, MPC_RNDNN);
   }
   else if (l == 59) {
      mpc_t z, tmp;
      mpc_init2 (z, mpc_get_prec (op));
      mpc_init2 (tmp, mpc_get_prec (rop));
      mpc_set (z, op, MPC_RNDNN);
      cm_modular_atkinhecke_eval (m, rop, z, 59, 5);
      cm_modular_atkinhecke_eval (m, tmp, z, 59, 29);
      mpc_add (rop, rop, tmp, MPC_RNDNN);
      mpc_add_ui (rop, rop, 1ul, MPC_RNDNN);
      mpc_clear (z);
      mpc_clear (tmp);
   }
   else if (l == 71) {
      mpc_t z, tmp;
      mpc_init2 (z, mpc_get_prec (op));
      mpc_init2 (tmp, mpc_get_prec (rop));
      mpc_set (z, op, MPC_RNDNN);
      cm_modular_atkinhecke_eval (m, rop, z, 71, 5);
      cm_modular_atkinhecke_eval (m, tmp, z, 71, 29);
      mpc_add (rop, rop, tmp, MPC_RNDNN);
      mpc_add_ui (rop, rop, 1ul, MPC_RNDNN);
      mpc_clear (z);
      mpc_clear (tmp);
   }
   else if (l == 131) {
      cm_modular_atkinhecke_eval (m, rop, op, 131, 61);
      mpc_add_ui (rop, rop, 1ul, MPC_RNDNN);
   }
   else {
      printf ("*** Called cm_modular_atkinhecke_level_eval with level %li, "
         "for which the optimal Atkin invariant is not implemented.\n", l);
      exit (1);
   }

}

/*****************************************************************************/
