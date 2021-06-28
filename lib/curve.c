/*

curve.c - code for computing cm curves

Copyright (C) 2009, 2010, 2021 Andreas Enge

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

static void elliptic_curve_double (mpz_ptr P_x, mpz_ptr P_y, bool *P_infty,
   mpz_srcptr a, mpz_srcptr p);
static void elliptic_curve_add (mpz_ptr P_x,  mpz_ptr P_y, bool *P_infty,
   mpz_srcptr Q_x, mpz_srcptr Q_y, const bool Q_infty,
   mpz_srcptr a, mpz_srcptr p);
static void elliptic_curve_multiply (mpz_ptr P_x, mpz_ptr P_y,
   bool *P_infty, mpz_srcptr m, mpz_srcptr a, mpz_srcptr p);
static void elliptic_curve_random (mpz_ptr P_x, mpz_ptr P_y,
   mpz_srcptr cofactor, mpz_srcptr a, mpz_srcptr b, mpz_srcptr p);
static bool curve_is_crypto (mpz_ptr l, mpz_ptr c, mpz_srcptr n,
   int_cl_t d, mpz_srcptr p, bool verbose);

/*****************************************************************************/

static void elliptic_curve_double (mpz_ptr P_x, mpz_ptr P_y, bool *P_infty,
   mpz_srcptr a, mpz_srcptr p)
   /* Replace P by 2P on the elliptic curve given by a over the prime field
      of characteristic p; b is implicit since it is assumed that the point
      lies on the curve.
      P_x and P_y are the X- and Y-coordinates of the point, respectively,
      P_infty is true if the point is the point at infinity. In this case,
      the values of P_x and P_y are not defined. */
{
   mpz_t factor, x_local, tmp, tmp2;

   mpz_init (factor);
   mpz_init (x_local);
   mpz_init (tmp);
   mpz_init (tmp2);

   if (!(*P_infty)) {
      /* There is nothing to do when P is already infinity */
      if (mpz_cmp_ui (P_y, 0ul) == 0)
         *P_infty = true;
      else {
         /* factor = (3*P_x^2 + a) / 2*P_y */
         mpz_pow_ui (tmp, P_x, 2ul);
         mpz_mod (tmp, tmp, p);
         mpz_mul_ui (factor, tmp, 3ul);
         mpz_add (factor, factor, a);
         mpz_mul_2exp (tmp, P_y, 1ul);
         mpz_invert (tmp2, tmp, p);
         mpz_mul (tmp, factor, tmp2);
         mpz_mod (factor, tmp, p);
         /* P_x = factor^2 - 2 P_x; save P_x */
         mpz_set (x_local, P_x);
         mpz_pow_ui (P_x, factor, 2ul);
         mpz_submul_ui (P_x, x_local, 2ul);
         mpz_mod (P_x, P_x, p);
         /* P_y = factor * (old x - new x) - P_y */
         mpz_sub (tmp, x_local, P_x);
         mpz_mul (tmp2, factor, tmp);
         mpz_sub (P_y, tmp2, P_y);
         mpz_mod (P_y, P_y, p);
      }
   }

   mpz_clear (factor);
   mpz_clear (x_local);
   mpz_clear (tmp);
   mpz_clear (tmp2);
}

/*****************************************************************************/

static void elliptic_curve_add (mpz_ptr P_x, mpz_ptr P_y, bool *P_infty,
   mpz_srcptr Q_x, mpz_srcptr Q_y, const bool Q_infty, mpz_srcptr a,
   mpz_srcptr p)
   /* Replace P by P+Q on the elliptic curve given by a and the implicit b
      over the prime field of characteristic p. */
{
   mpz_t factor, x_local, tmp, tmp2;

   mpz_init (factor);
   mpz_init (x_local);
   mpz_init (tmp);
   mpz_init (tmp2);

   if (!Q_infty) {
      if (*P_infty) {
         mpz_set (P_x, Q_x);
         mpz_set (P_y, Q_y);
         *P_infty = false;
      }
      else if (mpz_cmp (P_x, Q_x) == 0) {
         if (mpz_cmp (P_y, Q_y) == 0)
            elliptic_curve_double (P_x, P_y, P_infty, a, p);
         else
            *P_infty = true;
      }
      else {
         /* factor = Delta y / Delta x */
         mpz_sub (tmp, Q_x, P_x);
         mpz_invert (tmp2, tmp, p);
         mpz_sub (tmp, Q_y, P_y);
         mpz_mul (factor, tmp, tmp2);
         mpz_mod (factor, factor, p);
         /* P_x = factor^2 - (P_x + Q_x), keep a copy of P_x */
         mpz_set (x_local, P_x);
         mpz_pow_ui (tmp, factor, 2ul);
         mpz_add (tmp2, P_x, Q_x);
         mpz_sub (P_x, tmp, tmp2);
         mpz_mod (P_x, P_x, p);
         /* P_y = factor * (old P_x - new x) - P_y */
         mpz_sub (tmp, x_local, P_x);
         mpz_mul (tmp2, factor, tmp);
         mpz_sub (P_y, tmp2, P_y);
         mpz_mod (P_y, P_y, p);
      }
   }

   mpz_clear (factor);
   mpz_clear (x_local);
   mpz_clear (tmp);
   mpz_clear (tmp2);
}

/*****************************************************************************/

static void elliptic_curve_multiply (mpz_ptr P_x, mpz_ptr P_y, bool *P_infty,
   mpz_srcptr m, mpz_srcptr a, mpz_srcptr p)
   /* Replace P by mP on the elliptic curve given by a and the implicit b
      over the prime field of characteristic p. */
{
   mpz_t Q_x, Q_y, m_local, m_new;
   bool  Q_infty = true;
      /* m P is stored in Q; we use a right to left binary exponentiation    */
      /* scheme with the invariant Q + m P = const.                          */

   mpz_init (Q_x);
   mpz_init (Q_y);
   mpz_init (m_local);
   mpz_init (m_new);

   mpz_set (m_local, m);
   if (mpz_cmp_ui (m_local, 0ul) == 0)
      *P_infty = true;
   else if (!(*P_infty)) {
      if (mpz_cmp_ui (m_local, 0ul) < 0)
      {
         mpz_neg (m_local, m_local);
         mpz_neg (P_y, P_y);
      }
      while (mpz_cmp_ui (m_local, 0ul) > 0)
      {
         while (mpz_tdiv_q_ui (m_new, m_local, 2ul) == 0)
         {
            mpz_set (m_local, m_new);
            elliptic_curve_double (P_x, P_y, P_infty, a, p);
         }
         mpz_sub_ui (m_local, m_local, 1ul);
         elliptic_curve_add (Q_x, Q_y, &Q_infty, P_x, P_y, *P_infty, a, p);
      }
      /* copy the result */
      mpz_set (P_x, Q_x);
      mpz_set (P_y, Q_y);
      *P_infty = Q_infty;
   }

   mpz_clear (Q_x);
   mpz_clear (Q_y);
   mpz_clear (m_local);
   mpz_clear (m_new);
}

/*****************************************************************************/

static void elliptic_curve_random (mpz_ptr P_x, mpz_ptr P_y,
   mpz_srcptr cofactor, mpz_srcptr a, mpz_srcptr b, mpz_srcptr p)
   /* Create a point on the elliptic curve given by a and b over the prime
      field of characteristic p and multiply it by the cofactor until
      the result is different from infinity. If the curve order is cofactor
      times a prime, this results in a point of order this prime.
      The point is not really random, since successive X-coordinates from
      1 on are tested. */
{
   mpz_t  tmp;
   long unsigned int P_x_long = 0;
   bool P_infty = true;

   mpz_init (tmp);
   while (P_infty) {
      P_x_long++;
      /* P_y = P_x^3 + a P_x + b */
      mpz_mul_ui (P_y, a, P_x_long);
      mpz_add (P_y, P_y, b);
      mpz_add_ui (P_y, P_y, P_x_long * P_x_long * P_x_long);
      mpz_mod (P_y, P_y, p);
      /* try to compute the square root of P_y */
      if (mpz_jacobi (P_y, p) != -1) {
         mpz_set_ui (P_x, P_x_long);
         cm_nt_mpz_tonelli (P_y, P_y, p);
         /* get rid of the cofactor */
         P_infty = false;
         elliptic_curve_multiply (P_x, P_y, &P_infty, cofactor, a, p);
      }
   }
   mpz_clear (tmp);
}

/*****************************************************************************/

static bool curve_is_crypto (mpz_ptr l, mpz_ptr c, mpz_srcptr n,
   int_cl_t d, mpz_srcptr p, bool verbose)
   /* checks whether n might be a cryptographically secure cardinality for a */
   /* curve over F_p with discriminant d                                     */
   /* first tests if n, divided by a small cofactor, becomes a prime; if     */
   /* yes, the prime is returned via l, and the cofactor via c               */
   /* The cofactor which must be admitted depends on the discriminant;       */
   /* 4 : d = 1 (8); 4 | d and d/4 = 0 or 1 (4)                              */
   /* 2 : 4 | d and d/4 = 2 or 3 (4);                                        */
   /* 1 : d = 5 (8)                                                          */
   /* Watch out: Divisibility by the cofactor itself is not tested!          */
   /* Then tests for supersingular, anomalous and MOV curves.                */

{
   mpz_t t, rem;
   int k;

   if (d % 4 == 0)
   {
      if (d % 16 == 0 || ((d / 4) - 1) % 4 == 0)
      {
         mpz_tdiv_q_2exp (l, n, 2);
         mpz_set_ui (c, 4);
      }
      else
      {
         mpz_tdiv_q_2exp (l, n, 1);
         mpz_set_ui (c, 2);
      }
   }
   else if ((d - 1) % 8 == 0)
   {
      mpz_tdiv_q_2exp (l, n, 2);
      mpz_set_ui (c, 4);
   }
   else
   {
      mpz_set (l, n);
      mpz_set_ui (c, 1);
   }

   if (!cm_nt_is_prime (l)) {
      if (verbose)
         printf (".");
      return false;
   }

   mpz_init (t);
   mpz_init (rem);
   mpz_sub (t, n, p);
   mpz_sub_ui (t, t, 1);
   mpz_abs (t, t);

   if (mpz_cmp_ui (t, 0) == 0)
   {
      printf ("S");
      return false;
   }
   else if (mpz_cmp_ui (t, 1) == 0)
   {
      printf ("A");
      return false;
   }

   /* testing the MOV condition that q does not divide p^k - 1 for small */
   /* values of k                                                        */
   k = 1;
   mpz_set (t, p);
   mpz_mod (rem, t, l);
   while (k <= 9 && mpz_cmp_ui (rem, 1) != 0)
   {
      k++;
      mpz_mul (t, t, p);
      mpz_mod (rem, t, l);
   }
   if (k <= 9)
   {
      printf ("M%i", k);
      return false;
   }

   mpz_clear (t);
   mpz_clear (rem);

   return true;
}

/*****************************************************************************/

void cm_curve_crypto_param (mpz_ptr p, mpz_ptr n, mpz_ptr l, mpz_ptr c,
      int_cl_t d, int fieldsize, bool verbose)
   /* Given a discriminant d and a desired field size in bits (twice the
      bit security of the elliptic curve cryptosystem), return the
      cardinality of a cryptographically suitable elliptic curve.
      Precisely, p is the cardinality of the prime field, n the cardinality
      of the curve, l the prime order of a point on the curve and c=n/l the
      minimally possible cofactor for d, given as follows:
      d = 5 (8): u odd, v odd  (cofactor 1)
      d = 1 (8): u = 2 (4), v = 0 (4) (cofactor 4)
              or u = 0 (4), v = 2 (4) (cofactor 4, sometimes even higher)
      4 | d, d/4 = 2 (4): u = 2 (4), v odd (cofactor 2)
      4 | d, d/4 = 3 (4): u = 0 (4), v odd (cofactor 2)
      4 | d, d/4 = 1 (4): u = 2 (4), v even (cofactor 4)
                       or u = 0 (4), v odd (cofactor 8)
      4 | d, d/4 = 0 (4): u = 2 (4) (cofactor 4)
      To simplify the implementation, we step through the u and v in steps
      of 4. */
{
   mpz_t u, v, tmp;
   long unsigned int v_start;
   bool found = false;
   int deltav = 1000;

   if (fieldsize % 2 != 0)
      fieldsize++;

   mpz_init (u);
   mpz_init (v);
   mpz_init (tmp);

   mpz_set_ui (u, 1);
   mpz_mul_2exp (u, u, (fieldsize + 2) / 4);
   mpz_sub_ui (u, u, 2);
   mpz_mul_2exp (u, u, (fieldsize + 4) / 4);
   v_start = 1;

   /* so far, u is divisible by 4; update */
   if ((d - 5) % 8 == 0)
      mpz_add_ui (u, u, 1);
   else if (d % 8 == 0)
      mpz_add_ui (u, u, 2);
   else if ((d - 1) % 8 == 0 || (d % 4 == 0 && (d / 4 - 1) % 4 == 0))
   {
      mpz_add_ui (u, u, 2);
      v_start = 4;
   }

   while (!found)
   {
      if (deltav == 1000)
      {
         mpz_add_ui (u, u, 4);
         mpz_set_ui (v, v_start);
         deltav = 0;
      }
      else
      {
         mpz_add_ui (v, v, 4);
         deltav++;
      }

      mpz_mul (tmp, u, u);
      mpz_pow_ui (p, v, 2);
      mpz_mul_si (p, p, d);
      mpz_sub (p, tmp, p);
      /* should be divisible by 4... */
      mpz_tdiv_q_2exp (p, p, 2);
      if (cm_nt_is_prime (p))
      {
         mpz_add_ui (n, p, 1);
         mpz_sub (n, n, u);
         if (curve_is_crypto (l, c, n, d, p, verbose))
            found = true;
         else
         {
            mpz_add_ui (n, p, 1);
            mpz_add (n, n, u);
            if (curve_is_crypto (l, c, n, d, p, verbose))
               found = true;
         }
      }
   }

   if (verbose) {
      printf ("p   = "); mpz_out_str (stdout, 10, p); printf ("\n");
      printf ("u   = "); mpz_out_str (stdout, 10, u); printf ("\n");
      printf ("v   = "); mpz_out_str (stdout, 10, v); printf ("\n");
      printf ("n   = "); mpz_out_str (stdout, 10, n); printf ("\n");
      printf ("l   = "); mpz_out_str (stdout, 10, l); printf ("\n");
      printf ("N/l = "); mpz_out_str (stdout, 10, c); printf ("\n");
   }

   mpz_clear (u);
   mpz_clear (v);
   mpz_clear (tmp);
}

/*****************************************************************************/

void cm_curve_and_point (mpz_ptr a, mpz_ptr b, mpz_ptr x, mpz_ptr y,
   cm_param_srcptr param, cm_class_srcptr c,
   mpz_srcptr p, mpz_srcptr l, mpz_srcptr co,
   const char* modpoldir, bool verbose)
   /* Given CM parameters param, a class polynomial or class field tower
      stored in c, and curve cardinality parameters p (>=5, the cardinality
      of the prime field), a prime order l and a cofactor co, return curve
      parameters a and b defining an elliptic curve over F_p of cardinality
      n = l*co and a point P=(x,y) on the curve of order l.
      The algorithm will work in a slightly more general context
      (l and c are coprime, and gcd (exponent of curve group, l^\infty)=l),
      but the situation above is the common case for getting crypto curves
      or for ECPP. */
{
   mpz_t *j, *A, *B;
   mpz_t  twister;
   mpz_t  e, tmp, P_x, P_y;
   bool   P_infty;
   int    i, k, no_twists, no_j;
   bool   ok;
   cm_timer  clock;

   mpz_init (tmp);
   mpz_init (e);
   mpz_init (P_x);
   mpz_init (P_y);
   mpz_init (twister);

   /* Compute twister as a generator of F_p^* / (F_p^*)^n, where n is
      2, 4 or 6; this is equivalent to being a non-square such that
      additionally for n=6 (when p=1 mod 3), twister^((p-1)/3) != 1 mod p. */
   if (param->d != -3) {
      mpz_set_ui (twister, 2);
      while (mpz_jacobi (twister, p) != -1)
         mpz_add_ui (twister, twister, 1);
      no_twists = (param->d == -4 ? 4 : 2);
   }
   else {
      mpz_sub_ui (e, p, 1);
      if (!mpz_divisible_ui_p (e, 3)) {
         printf ("*** Error: p != 1 mod 3 for d=-3\n");
         exit (1);
      }
      mpz_divexact_ui (e, e, 3);
      mpz_set_ui (twister, 1);
      do {
         mpz_add_ui (twister, twister, 1);
         mpz_powm (tmp, twister, e, p);
      } while (mpz_jacobi (twister, p) != -1 || !mpz_cmp_ui (tmp, 1));
      no_twists = 6;
   }

   A = (mpz_t *) malloc (no_twists * sizeof (mpz_t));
   B = (mpz_t *) malloc (no_twists * sizeof (mpz_t));
   for (i = 0; i < no_twists; i++) {
      mpz_init (A [i]);
      mpz_init (B [i]);
   }

   cm_timer_continue (cm_timer1);
   j = cm_class_get_j_mod_p (&no_j, param, c, p, modpoldir, verbose);
   cm_timer_stop (cm_timer1);

   cm_timer_start (clock);
   cm_timer_continue (cm_timer2);
   ok = false;
   for (i = 0; i < no_j && !ok; i++) {
      /* Construct one curve with the given j-invariant. */
      if (mpz_cmp_ui (j [i], 1728) == 0) {
         mpz_set_ui (A [0], 1);
         mpz_set_ui (B [0], 0);
      }
      else if (mpz_cmp_ui (j [i], 0) == 0) {
         mpz_set_ui (A [0], 0);
         mpz_set_ui (B [0], 1);
      }
      else {
         /* a = 3 * j [i] * (1728 - j [i]),
            b = 2 * j [i] * (1728 - j [i])^2 */
         mpz_ui_sub (e, 1728, j [i]);
         mpz_mul (A [0], j [i], e);
         mpz_mod (A [0], A [0], p);
         mpz_mul (B [0], A [0], e);
         mpz_mul_ui (A [0], A [0], 3);
         mpz_mod (A [0], A [0], p);
         mpz_mul_2exp (B [0], B [0], 1);
         mpz_mod (B [0], B [0], p);
      }

      /* Compute the twists. */
      if (no_twists == 2) {
         mpz_powm_ui (tmp, twister, 2, p);
         mpz_mul (A [1], A [0], tmp);
         mpz_mod (A [1], A [1], p);
         mpz_mul (tmp, tmp, twister);
         mpz_mod (tmp, tmp, p);
         mpz_mul (B [1], B [0], tmp);
         mpz_mod (B [1], B [1], p);
      }
      else if (no_twists == 4)
         for (k = 1; k < no_twists; k++) {
            mpz_mul (A [k], A [k-1], twister);
            mpz_mod (A [k], A [k], p);
            mpz_set_ui (B [k], 0);
         }
      else
         for (k = 1; k < no_twists; k++) {
            mpz_mul (B [k], B [k-1], twister);
            mpz_mod (B [k], B [k], p);
            mpz_set_ui (A [k], 0);
         }

      /* Go through the twists and look for a suitable curve. */
      for (k = 0; k < no_twists && !ok; k++) {
         mpz_set (a, A [k]);
         mpz_set (b, B [k]);
         cm_timer_continue (cm_timer3);
         elliptic_curve_random (P_x, P_y, co, a, b, p);
         cm_timer_stop (cm_timer3);
         mpz_set (x, P_x);
         mpz_set (y, P_y);
         P_infty = false;
         cm_timer_continue (cm_timer4);
         elliptic_curve_multiply (P_x, P_y, &P_infty, l, a, p);
         cm_timer_stop (cm_timer4);
         if (P_infty)
            ok = true;
      }
   }
   cm_timer_stop (cm_timer2);

   if (!ok) {
      printf ("\n*** No suitable curve found!\n");
      exit (1);
   }

   cm_timer_stop (clock);
   if (verbose) {
      printf ("p = "); mpz_out_str (stdout, 10, p); printf ("\n");
      printf ("n = "); mpz_out_str (stdout, 10, co);
      printf (" * "); mpz_out_str (stdout, 10, l); printf ("\n");
      printf ("j = "); mpz_out_str (stdout, 10, j [i-1]); printf ("\n");
      printf ("a = "); mpz_out_str (stdout, 10, a); printf ("\n");
      printf ("b = "); mpz_out_str (stdout, 10, b); printf ("\n");
      printf ("--- Time for curve: %.1f\n", cm_timer_get (clock));
   }

   for (i = 0; i < no_j; i++)
      mpz_clear (j [i]);
   free (j);
   for (i = 0; i < no_twists; i++) {
      mpz_clear (A [i]);
      mpz_clear (B [i]);
   }
   free (A);
   free (B);
   mpz_clear (tmp);
   mpz_clear (e);
   mpz_clear (P_x);
   mpz_clear (P_y);
   mpz_clear (twister);
}

/*****************************************************************************/
/*****************************************************************************/
