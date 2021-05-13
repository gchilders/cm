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

static void elliptic_curve_double (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t a, mpz_t p);
static void elliptic_curve_add (mpz_t P_x,  mpz_t P_y, bool *P_infty,
   mpz_t Q_x, mpz_t Q_y, const bool Q_infty, mpz_t a, mpz_t p);
static void elliptic_curve_multiply (mpz_t P_x, mpz_t P_y,
   bool *P_infty, mpz_t m, mpz_t a, mpz_t p);
static void elliptic_curve_random (mpz_t P_x, mpz_t P_y,
   mpz_t cofactor, mpz_t a, mpz_t b, mpz_t p);
static bool curve_is_crypto (mpz_t l, mpz_t c, mpz_t n, int_cl_t d,
   mpz_t p, bool verbose);
static void curve_compute_param (mpz_t p, mpz_t n, mpz_t l, mpz_t c,
      int_cl_t d, int fieldsize, bool verbose);

/*****************************************************************************/

static void elliptic_curve_double (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t a, mpz_t p)
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

static void elliptic_curve_add (mpz_t P_x,  mpz_t P_y, bool *P_infty,
   mpz_t Q_x, mpz_t Q_y, const bool Q_infty, mpz_t a, mpz_t p)
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

static void elliptic_curve_multiply (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t m, mpz_t a, mpz_t p)
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
   *P_infty = false;

   mpz_set (m_local, m);
   if (mpz_cmp_ui (m_local, 0ul) == 0)
      *P_infty = true;
   else
   {
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

static void elliptic_curve_random (mpz_t P_x, mpz_t P_y,
   mpz_t cofactor, mpz_t a, mpz_t b, mpz_t p)
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
         cm_nt_mpz_tonelli_z (P_y, P_y, p);
         /* get rid of the cofactor */
         P_infty = false;
         elliptic_curve_multiply (P_x, P_y, &P_infty, cofactor, a, p);
      }
   }
   mpz_clear (tmp);
}

/*****************************************************************************/

static bool curve_is_crypto (mpz_t l, mpz_t c, mpz_t n, int_cl_t d,
   mpz_t p, bool verbose)
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

static void curve_compute_param (mpz_t p, mpz_t n, mpz_t l, mpz_t c,
      int_cl_t d, int fieldsize, bool verbose)
   /* computes the curve parameters p (size of prime field), n (cardinality  */
   /* of curve), l (prime order of point) and c (n/l)                        */
   /* Try u and v until u^2 - v^2 d is 4 times a prime. Congruences of u and */
   /* v modulo 2 and 4 are taken into account, also to make sure that the    */
   /* cofactor of the curve cardinality can be as small as possible.         */
   /* Conditions:                                                            */
   /* d = 5 (8): u odd, v odd  (cofactor 1)                                  */
   /* d = 1 (8): u = 2 (4), v = 0 (4) (cofactor 4)                           */
   /*         or u = 0 (4), v = 2 (4) (cofactor 4, sometimes even higher)    */
   /* 4 | d, d/4 = 2 (4): u = 2 (4), v odd (cofactor 2)                      */
   /* 4 | d, d/4 = 3 (4): u = 0 (4), v odd (cofactor 2)                      */
   /* 4 | d, d/4 = 1 (4): u = 2 (4), v even (cofactor 4)                     */
   /*                  or u = 0 (4), v odd (cofactor 8)                      */
   /* 4 | d, d/4 = 0 (4): u = 2 (4) (cofactor 4)                             */
   /* To simplify the implementation, we step through the u and v in steps   */
   /* of 4.                                                                  */
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

void cm_curve_compute_curve (int_cl_t d, char inv, int fieldsize,
   const char* modpoldir, bool pari, bool readwrite, bool print,
   bool tower, bool verbose)
   /* computes a curve with the good number of points, using the given       */
   /* invariant                                                              */
   /* pari indicates whether the pari stack should be initialised; if set
      to false, this is the responsability of the calling function.
      If readwrite is true, tries to read the polynomial from a file, and    */
   /* computes it and writes it to the file if it is not yet present.        */
   /* print indicates whether the curve parameters shall be displayed on     */
   /* screen                                                                 */
   /* If tower is set, go through a class field tower instead of a class
      polynomial. */
{
   mpz_t  p;
      /* the size of the base field */
   mpz_t  n, l, c;
      /* the cardinality of the curve, the prime order of a base point and  */
      /* the cofactor C = N/Q                                               */
   mpz_t* j;
   mpz_t  nonsquare;
      /* a quadratic non-residue modulo P to compute quadratic twists        */
   mpz_t  a, b, k, P_x, P_y, x, y;
      /* the curve parameters mod P, an auxiliary variable and point         */
      /* coordinates of a random point on the curve (twice)                  */
   mpz_t  tmp;
   bool   P_infty;
   int    i, no_j;
   bool   ok = false;
      /* true if at least one suitable curve could be found */
   int    count = 0;
      /* number of suitable curves */
   cm_timer  clock;

   mpz_init (p);
   mpz_init (n);
   mpz_init (l);
   mpz_init (c);
   mpz_init (tmp);
   mpz_init (a);
   mpz_init (b);
   mpz_init (k);
   mpz_init (P_x);
   mpz_init (P_y);
   mpz_init (x);
   mpz_init (y);
   mpz_init_set_ui (nonsquare, 2);

   cm_timer_start (clock);
   curve_compute_param (p, n, l, c, d, fieldsize, verbose);
   cm_timer_stop (clock);
   while (mpz_jacobi (nonsquare, p) != -1)
      mpz_add_ui (nonsquare, nonsquare, 1);
   if (verbose)
      printf ("--- Time for P: %.1f\n\n", cm_timer_get (clock));

   j = cm_class_get_j_mod_P (d, inv, p, &no_j, modpoldir, pari, readwrite,
      tower, verbose);

   cm_timer_start (clock);
   for (i = 0; i < no_j && !ok; i++)
   {
      /* construct one curve with the given j-invariant */
      if (mpz_cmp_ui (j [i], 1728) == 0)
      /* may occur with spurious roots of modular polynomials */
      {
         mpz_set_ui (a, 1);
         mpz_set_ui (b, 0);
      }
      else if (mpz_cmp_ui (j [i], 0) == 0)
      /* should not occur as d \not= -3 */
      {
         mpz_set_ui (a, 0);
         mpz_set_ui (b, 1);
      }
      else
      {
         /* k = j [i] / (1728 - j_mod) */
         mpz_sub_ui (k, j [i], 1728);
         mpz_neg (k, k);
         mpz_invert (tmp, k, p);
         mpz_mul (k, tmp, j [i]);
         /* b = 2 * k; */
         mpz_mul_2exp (b, k, 1);
         mpz_mod (b, b, p);
         /* a = 3 * k; */
         mpz_add (a, b, k);
         mpz_mod (a, a, p);
      }

      elliptic_curve_random (P_x, P_y, c, a, b, p);
      mpz_set (x, P_x);
      mpz_set (y, P_y);
      P_infty = false;
      elliptic_curve_multiply (P_x, P_y, &P_infty, l, a, p);
      if (P_infty) {
         ok = true;
         count++;
      }
      else {
         /* consider the quadratic twist */
         mpz_pow_ui (tmp, nonsquare, 2);
         mpz_mul (a, a, tmp);
         mpz_mod (a, a, p);
         mpz_mul (b, b, tmp);
         mpz_mul (b, b, nonsquare);
         mpz_mod (b, b, p);
         elliptic_curve_random (P_x, P_y, c, a, b, p);
         mpz_set (x, P_x);
         mpz_set (y, P_y);
         P_infty = false;
         elliptic_curve_multiply (P_x, P_y, &P_infty, l, a, p);
         if (P_infty) {
            ok = true;
            count++;
         }
      }
   }

   if (!ok) {
      printf ("\n*** No suitable curve found!\n");
      exit (1);
   }

   cm_timer_stop (clock);
   if (print) {
      printf ("p = "); mpz_out_str (stdout, 10, p); printf ("\n");
      printf ("n = "); mpz_out_str (stdout, 10, n); printf ("\n");
      printf ("  = "); mpz_out_str (stdout, 10, c);
      printf (" * "); mpz_out_str (stdout, 10, l); printf ("\n");
      printf ("j = "); mpz_out_str (stdout, 10, j [i-1]); printf ("\n");
      printf ("a = "); mpz_out_str (stdout, 10, a); printf ("\n");
      printf ("b = "); mpz_out_str (stdout, 10, b); printf ("\n");
   }
   if (verbose)
      printf ("--- Time for curve: %.1f\n", cm_timer_get (clock));

   for (i = 0; i < no_j; i++)
      mpz_clear (j [i]);
   free (j);
   mpz_clear (p);
   mpz_clear (n);
   mpz_clear (l);
   mpz_clear (c);
   mpz_clear (tmp);
   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (k);
   mpz_clear (P_x);
   mpz_clear (P_y);
   mpz_clear (x);
   mpz_clear (y);
   mpz_clear (nonsquare);
}

/*****************************************************************************/
/*****************************************************************************/
