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

void gw_powm ( mpz_ptr mr, mpz_ptr mb_inp, mpz_ptr me, mpz_srcptr mn );
static bool elliptic_curve_dehomogenise (mpz_ptr x, mpz_ptr y,
   mpz_ptr z, mpz_srcptr p);
static void elliptic_curve_double (mpz_ptr x, mpz_ptr y, mpz_ptr z,
   mpz_ptr t, mpz_srcptr p);
static void elliptic_curve_mixedadd (mpz_ptr x1, mpz_ptr y1, mpz_ptr z1,
   mpz_ptr t1, mpz_srcptr x2, mpz_srcptr y2, mpz_srcptr a, mpz_srcptr p);
static void elliptic_curve_multiply (mpz_ptr P_x, mpz_ptr P_y,
   bool *P_infty, mpz_srcptr m, mpz_srcptr a, mpz_srcptr p);
static void elliptic_curve_random (mpz_ptr P_x, mpz_ptr P_y,
   mpz_srcptr cofactor, mpz_srcptr a, mpz_srcptr b, mpz_srcptr p);
static bool curve_is_crypto (mpz_ptr l, mpz_ptr c, mpz_srcptr n,
   int_cl_t d, mpz_srcptr p, bool verbose);

/*****************************************************************************/

static bool elliptic_curve_dehomogenise (mpz_ptr x, mpz_ptr y,
   mpz_ptr z, mpz_srcptr p)
   /* If the point in Jacobian coordinates is infinite, return true;
      otherwise replace it by the corresponding affine point with z=1
      and return false. */
{
   mpz_t tmp;

   if (!mpz_cmp_ui (z, 0ul))
      return true;
   else {
      mpz_init (tmp);

      mpz_invert (z, z, p);
      mpz_mul (tmp, z, z);
      mpz_mod (tmp, tmp, p);
      mpz_mul (x, x, tmp);
      mpz_mod (x, x, p);
      mpz_mul (tmp, tmp, z);
      mpz_mod (tmp, tmp, p);
      mpz_mul (y, y, tmp);
      mpz_mod (y, y, p);
      mpz_set_ui (z, 1ul);

      mpz_clear (tmp);

      return false;
   }
}
 
/*****************************************************************************/

static void elliptic_curve_double (mpz_ptr x, mpz_ptr y, mpz_ptr z,
   mpz_ptr t, mpz_srcptr p)
   /* Replace P, given in modified Jacobian coordinates (x:y:z:t), by 2P.
      The parameters a and b are implicit: a=z^4/t, and b since the
      point lies on the curve. */
{
   mpz_t xx, yy, yyyy, s, m, x2, y2, z2, t2, tmp;

   mpz_init (xx);
   mpz_init (yy);
   mpz_init (yyyy);
   mpz_init (s);
   mpz_init (m);
   mpz_init (x2);
   mpz_init (y2);
   mpz_init (z2);
   mpz_init (t2);
   mpz_init (tmp);

   if (mpz_cmp_ui (z, 0ul)) {
      /* There is nothing to do when P is already infinity.
         Otherwise we use the formula of
         https://www.hyperelliptic.org/EFD/g1p/auto-shortw-modified.html#doubling-dbl-2009-bl
         that requires 3 multiplications and 5 squarings; we also try to
         minimise the reductions modulo p, carrying them out only when an
         element is used in a subsequent multiplication or is a part of the
         result, and use 8 reductions. */
      /* xx = x^2, yy = 2*y^2, yyyy = 4*y^4 */
      mpz_mul (xx, x, x);
      mpz_mul (yy, y, y);
      mpz_mul_2exp (yy, yy, 1);
      mpz_mod (yy, yy, p);
      mpz_mul (yyyy, yy, yy);
      mpz_mod (yyyy, yyyy, p);
      /* s = (x + yy)^2 - xx - yyyy = 8*x*y^2 */
      mpz_add (s, x, yy);
      mpz_mul (s, s, s);
      mpz_sub (s, s, xx);
      mpz_sub (s, s, yyyy);
      mpz_mod (s, s, p);
      /* m = 3 * xx + t */
      mpz_mul_ui (m, xx, 3);
      mpz_add (m, m, t);
      mpz_mod (m, m, p);
      /* x2 = m^2 - 2 * s */
      mpz_mul (x2, m, m);
      mpz_mul_2exp (tmp, s, 1);
      mpz_sub (x2, x2, tmp);
      mpz_mod (x2, x2, p);
      /* y2 = m * (s - x2) - 2 * yyyy */
      mpz_sub (y2, s, x2);
      mpz_mul (y2, y2, m);
      mpz_mul_2exp (tmp, yyyy, 1);
      mpz_sub (y2, y2, tmp);
      mpz_mod (y2, y2, p);
      /* z2 = 2*y*z */
      mpz_mul (z2, y, z);
      mpz_mul_2exp (z2, z2, 1);
      mpz_mod (z2, z2, p);
      /* t2 = 4*yyyy*t = 16*y^4*t */
      mpz_mul (t2, yyyy, t);
      mpz_mul_2exp (t2, t2, 2);
      mpz_mod (t2, t2, p);

      mpz_set (x, x2);
      mpz_set (y, y2);
      mpz_set (z, z2);
      mpz_set (t, t2);
   }

   mpz_clear (xx);
   mpz_clear (yy);
   mpz_clear (yyyy);
   mpz_clear (s);
   mpz_clear (m);
   mpz_clear (x2);
   mpz_clear (y2);
   mpz_clear (z2);
   mpz_clear (t2);
   mpz_clear (tmp);
}

/*****************************************************************************/

static void elliptic_curve_mixedadd (mpz_ptr x1, mpz_ptr y1, mpz_ptr z1,
   mpz_ptr t1, mpz_srcptr x2, mpz_srcptr y2, mpz_srcptr a, mpz_srcptr p)
   /* Replace P1, given by modified Jacobian coordinates (x1:y1:z1:t1),
      by P1+P2, where P2 is non-zero and given by affine coordinates
      (x2:y2:1:a), on the same elliptic curve given by a over the prime
      field of characteristic p. */
{
   mpz_t z1z1, h, hh, i, j, r, v, x3, y3, z3, t3, tmp;

   if (!mpz_cmp_ui (z1, 0ul)) {
      /* P1 is zero */
      mpz_set (x1, x2);
      mpz_set (y1, y2);
      mpz_set_ui (z1, 1ul);
      mpz_set (t1, a);
   }
   else {
      mpz_init (z1z1);
      mpz_init (h);
      mpz_init (r);

      /* Use the formula for mixed addition of
         https://www.hyperelliptic.org/EFD/g1p/auto-shortw-modified.html#addition-madd-2009-bl
         which requires 8 multiplications and 6 squarings, as well as
         13 reductions modulo p.
         On the way, we need to check whether P1 == +- P2, for which
         we slightly reorder the formula to see this happen earlier. */
      /* z1z1 = z1^2 */
      mpz_mul (z1z1, z1, z1);
      mpz_mod (z1z1, z1z1, p);
      /* h = x2*z1z1 - x1 */
      mpz_mul (h, x2, z1z1);
      mpz_sub (h, h, x1);
      mpz_mod (h, h, p);
      /* r = 2 * (y2 * z1^3 - y1) */
      mpz_mul (r, z1z1, z1);
      mpz_mod (r, r, p);
      mpz_mul (r, r, y2);
      mpz_sub (r, r, y1);
      mpz_mul_2exp (r, r, 1ul);
      mpz_mod (r, r, p);
      if (!mpz_cmp_ui (h, 0ul))
         /* P1 == +- P2 */
         if (!mpz_cmp_ui (r, 0ul))
            /* P1 == P2 */
            elliptic_curve_double (x1, y1, z1, t1, p);
         else
            /* P1 == -P2 */
            mpz_set_ui (z1, 0ul);
      else {
         mpz_init (hh);
         mpz_init (i);
         mpz_init (j);
         mpz_init (v);
         mpz_init (x3);
         mpz_init (y3);
         mpz_init (z3);
         mpz_init (t3);
         mpz_init (tmp);

         /* hh = h^2 */
         mpz_mul (hh, h, h);
         /* i = 4 * hh */
         mpz_mul_2exp (i, hh, 2);
         mpz_mod (i, i, p);
         /* j = h * i */
         mpz_mul (j, h, i);
         mpz_mod (j, j, p);
         /* v = x1 * i */
         mpz_mul (v, x1, i);
         mpz_mod (v, v, p);
         /* x3 = r^2 - j - 2 * v */
         mpz_mul (x3, r, r);
         mpz_sub (x3, x3, j);
         mpz_mul_2exp (tmp, v, 1);
         mpz_sub (x3, x3, tmp);
         mpz_mod (x3, x3, p);
         /* y3 = r * (v - x3) - 2 * y1 * j */
         mpz_sub (y3, v, x3);
         mpz_mul (y3, y3, r);
         mpz_mul (tmp, y1, j);
         mpz_mul_2exp (tmp, tmp, 1);
         mpz_sub (y3, y3, tmp);
         mpz_mod (y3, y3, p);
         /* z3 = (z1 + h)^2 - z1z1 - hh = 2 * z1 * h */
         mpz_add (z3, z1, h);
         mpz_mul (z3, z3, z3);
         mpz_sub (z3, z3, z1z1);
         mpz_sub (z3, z3, hh);
         mpz_mod (z3, z3, p);
         /* t3 = a*z3^4 */
         mpz_mul (t3, z3, z3);
         mpz_mod (t3, t3, p);
         mpz_mul (t3, t3, t3);
         mpz_mod (t3, t3, p);
         mpz_mul (t3, t3, a);
         mpz_mod (t3, t3, p);

         mpz_set (x1, x3);
         mpz_set (y1, y3);
         mpz_set (z1, z3);
         mpz_set (t1, t3);

         mpz_clear (hh);
         mpz_clear (i);
         mpz_clear (j);
         mpz_clear (v);
         mpz_clear (x3);
         mpz_clear (y3);
         mpz_clear (z3);
         mpz_clear (t3);
         mpz_clear (tmp);
      }

      mpz_clear (z1z1);
      mpz_clear (h);
      mpz_clear (r);
   }
}

/*****************************************************************************/

static void elliptic_curve_multiply (mpz_ptr P_x, mpz_ptr P_y, bool *P_infty,
   mpz_srcptr m, mpz_srcptr a, mpz_srcptr p)
   /* Replace P by mP on the elliptic curve given by a and the implicit b
      over the prime field of characteristic p, where m is assumed to be
      non-negative. */
{
   mpz_t x, y, z, t;
      /* m P is stored in modified Jacobian coordinates as (x:y:z:t);
         we use a left to right sliding window exponentiation scheme. */
   int w; /* window width */
   int l; /* length of precomputed table 2^(w-1) */
   mpz_t *mx, *my, *mz, *mt;
      /* (mx [i] : my [i] : mz [i] : mt [i]) contains the modified Jacobian
         coordinates of the multiple (2*i+1)*P for 0 <= i < l. These
         multiples are eventually dehomogenised, but care must be taken
         if by chance any of them is 0. */
   int i, n, j, start, stop;
   bool infty;

   if (!mpz_cmp_ui (m, 0ul))
      *P_infty = true;
   else if (!(*P_infty || !mpz_cmp_ui (m, 1ul))) {
      /* From now on P is a finite point and m at least 2.
         Precompute 2*P and the odd multiples of P, and handle the special
         cases when infinity occurs as one of them. */

      mpz_init_set (x, P_x);
      mpz_init_set (y, P_y);
      mpz_init_set_ui (z, 1ul);
      mpz_init_set (t, a);
      elliptic_curve_double (x, y, z, t, p);
      if (elliptic_curve_dehomogenise (x, y, z, p))
         /* P is of order 2. So we return either P or infinity depending
            on the parity of m. */
         *P_infty = mpz_divisible_2exp_p (m, 1);
      else {
         /* Now (x, y) == 2*P. Precompute the odd multiples of P. */

         /* For the window width, we use the following simplified model:
            For a width of k, the precomputations take 2^(k-1) additions;
            the average window length including the following run of 0 is
            k+1, so the total number of additions is
            f (n, k) = 2^(k-1) + n/(k+1) with n = ld(m).
            This function is unimodular, so we should set w to the smallest
            k such that f(n,k) <= f(n,k+1) <=> (k+1)*(k+2)*2^(k-1) <= n. */
         j = mpz_sizeinbase (m, 2);
         if (j <= 6)
            w = 1;
         else if (j <= 24)
            w = 2;
         else if (j <= 80)
            w = 3;
         else if (j <= 240)
            w = 4;
         else if (j <= 672)
            w = 5;
         else if (j <= 1792)
            w = 6;
         else if (j <= 4608)
            w = 7;
         else if (j <= 11520)
            w = 8;
         else if (j <= 28160)
            w = 9;
         else if (j <= 67584)
            w = 10;
         else
            w = 11;
         l = 1 << (w - 1);
         mx = (mpz_t *) malloc (l * sizeof (mpz_t));
         my = (mpz_t *) malloc (l * sizeof (mpz_t));
         mz = (mpz_t *) malloc (l * sizeof (mpz_t));
         mt = (mpz_t *) malloc (l * sizeof (mpz_t));
         for (i = 0; i < l; i++) {
            mpz_init (mx [i]);
            mpz_init (my [i]);
            mpz_init (mz [i]);
            mpz_init (mt [i]);
         }
         mpz_set (mx [0], P_x);
         mpz_set (my [0], P_y);
         mpz_set_ui (mz [0], 1ul);
         mpz_set (mt [0], a);
         infty = false;
         for (i = 1; !infty && i < l; i++) {
            mpz_set (mx [i], mx [i-1]);
            mpz_set (my [i], my [i-1]);
            mpz_set (mz [i], mz [i-1]);
            mpz_set (mt [i], mt [i-1]);
            elliptic_curve_mixedadd (mx [i], my [i], mz [i], mt [i],
               x, y, a, p);
            infty = !mpz_cmp_ui (mz [i], 0ul);
         }
         if (infty) {
            /* P has order 2*i-1.
               Let n = m % (2*i-1) = 2^j * (2*k + 1).
               Then m*P = 2^j * (mx [k] : my [k] : mz [k] : mt [k]) can be
               computed with j doublings. This is a special case of the
               sliding window code below with only one window. */
            n = mpz_fdiv_ui (m, (unsigned long int) (2*i-1));
            j = 0;
            while (n % 2 == 0) {
               n /= 2;
               j++;
            }
            n /= 2;
            mpz_set (x, mx [n]);
            mpz_set (y, my [n]);
            mpz_set (z, mz [n]);
            mpz_set (t, mt [n]);
            for (i = 0; i < j; i++)
               elliptic_curve_double (x, y, z, t, p);
         }
         else {
            /* None of the precomputed points is infinite, use
               sliding windows. */

            for (i = 1; i < l; i++)
               /* This could be done more efficiently by inverting all
                  Z-coordinates at the same time using Montgomery's trick,
                  but right now this is not a bottleneck. */
               elliptic_curve_dehomogenise (mx [i], my [i], mz [i], p);

            start = mpz_sizeinbase (m, 2) - 2;
               /* position of the next bit to treat */
            if (!mpz_tstbit (m, start))
               /* Microoptimisation: Start with 2*P = (x, y). */
               start--;
            else {
               mpz_set (x, P_x);
               mpz_set (y, P_y);
            }
            mpz_set_ui (z, 1ul);
            mpz_set (t, a);

            while (start >= 0) {
               /* Slide to the right while the bit is 0. */
               while (start >= 0 && !mpz_tstbit (m, start)) {
                  elliptic_curve_double (x, y, z, t, p);
                  start--;
               }
               if (start >= 0) {
                  /* Start is on a bit 1, find the end of the window on
                     another bit 1. */
                  stop = (start < w-1 ? 0 : start - (w - 1));
                  stop = mpz_scan1 (m, stop);
                  /* Double according to the length of the window. */
                  j = start - stop + 1;
                  for (i = 0; i < j; i++)
                     elliptic_curve_double (x, y, z, t, p);
                  /* Add the content of the window. */
                  if (stop == start)
                     j = 0;
                  else
                     for (j = 1, i = start - 1; i > stop; i--)
                        j = 2 * j + mpz_tstbit (m, i);
                  elliptic_curve_mixedadd (x, y, z, t, mx [j], my [j],
                     a, p);
                  start = stop - 1;
               }
            }

         }

         *P_infty = elliptic_curve_dehomogenise (x, y, z, p);
         mpz_set (P_x, x);
         mpz_set (P_y, y);

         for (i = 0; i < l; i++) {
            mpz_clear (mx [i]);
            mpz_clear (my [i]);
            mpz_clear (mz [i]);
            mpz_clear (mt [i]);
         }
         free (mx);
         free (my);
         free (mz);
         free (mt);
      }

      mpz_clear (x);
      mpz_clear (y);
      mpz_clear (z);
      mpz_clear (t);
   }
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

void cm_curve_and_point_stat (mpz_ptr a, mpz_ptr b, mpz_ptr x, mpz_ptr y,
   cm_param_srcptr param, cm_class_srcptr c,
   mpz_srcptr p, mpz_srcptr l, mpz_srcptr co,
   const char* modpoldir, bool print, bool verbose, cm_stat_t stat)
   /* Given CM parameters param, a class polynomial or class field tower
      stored in c, and curve cardinality parameters p (>=5, the cardinality
      of the prime field), a prime order l and a cofactor co, return curve
      parameters a and b defining an elliptic curve over F_p of cardinality
      n = l*co and a point P=(x,y) on the curve of order l.
      The algorithm will work in a slightly more general context
      (l and c are coprime, and gcd (exponent of curve group, l^\infty)=l),
      but the situation above is the common case for getting crypto curves
      or for ECPP.
      The parameter print indicates whether the resulting curve and point
      parameters are output on screen; verbose indicates whether additional
      information is printed during the execution. */
{
   mpz_t *j, *A, *B;
   mpz_t  twister;
   mpz_t  e, tmp, P_x, P_y;
   bool   P_infty;
   int    i, k, no_twists, no_j;
   bool   ok;
   cm_timer_t clock;

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
         gw_powm (tmp, twister, e, p);
      } while (mpz_jacobi (twister, p) != -1 || !mpz_cmp_ui (tmp, 1));
      no_twists = 6;
   }

   A = (mpz_t *) malloc (no_twists * sizeof (mpz_t));
   B = (mpz_t *) malloc (no_twists * sizeof (mpz_t));
   for (i = 0; i < no_twists; i++) {
      mpz_init (A [i]);
      mpz_init (B [i]);
   }

   if (stat != NULL)
      cm_timer_continue (stat->timer [2]);
   j = cm_class_get_j_mod_p (&no_j, param, c, p, modpoldir, verbose);
   if (stat != NULL)
      cm_timer_stop (stat->timer [2]);

   cm_timer_start (clock);
   if (stat !=NULL)
      cm_timer_continue (stat->timer [3]);
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
         elliptic_curve_random (P_x, P_y, co, a, b, p);
         mpz_set (x, P_x);
         mpz_set (y, P_y);
         P_infty = false;
         elliptic_curve_multiply (P_x, P_y, &P_infty, l, a, p);
         if (P_infty)
            ok = true;
      }
   }
   if (stat != NULL)
      cm_timer_stop (stat->timer [3]);

   if (!ok) {
      printf ("\n*** No suitable curve found!\n");
      exit (1);
   }

   cm_timer_stop (clock);
   if (print || verbose) {
      printf ("p = "); mpz_out_str (stdout, 10, p); printf ("\n");
      printf ("n = "); mpz_out_str (stdout, 10, co);
      printf (" * "); mpz_out_str (stdout, 10, l); printf ("\n");
      printf ("a = "); mpz_out_str (stdout, 10, a); printf ("\n");
      printf ("b = "); mpz_out_str (stdout, 10, b); printf ("\n");
      printf ("x = "); mpz_out_str (stdout, 10, x); printf ("\n");
      printf ("y = "); mpz_out_str (stdout, 10, y); printf ("\n");
   }
   if (verbose) {
      printf ("j = "); mpz_out_str (stdout, 10, j [i-1]); printf ("\n");
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

void cm_curve_and_point (mpz_ptr a, mpz_ptr b, mpz_ptr x, mpz_ptr y,
   cm_param_srcptr param, cm_class_srcptr c,
   mpz_srcptr p, mpz_srcptr l, mpz_srcptr co,
   const char* modpoldir, bool print, bool verbose)
{
   cm_curve_and_point_stat (a, b, x, y, param, c, p, l, co, modpoldir,
      print, verbose, NULL);
}

/*****************************************************************************/
/*****************************************************************************/
