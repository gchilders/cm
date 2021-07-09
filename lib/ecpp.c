/*

ecpp.c - code for computing ECPP certificates

Copyright (C) 2021 Andreas Enge

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

static void compute_h (uint_cl_t *h, uint_cl_t Dmax);
static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   int no_old, int no_new);
static int_cl_t** compute_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar, uint_cl_t Dmax, uint_cl_t hmaxprime,
   uint_cl_t *h);
static int disc_cmp (const void* d1, const void* d2);
static void trial_div (mpz_ptr l, mpz_srcptr n, mpz_srcptr primorialB);
static bool is_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   mpz_srcptr root, const unsigned int delta, mpz_srcptr primorialB,
   int_cl_t d);
static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const unsigned int delta, mpz_srcptr primorialB);

/*****************************************************************************/

static void compute_h (uint_cl_t *h, uint_cl_t Dmax)
   /* Assuming that h is an array of length (Dmax/2)+1 for even Dmax,
      compute in h [|D|/2] the class number for fundamental discriminants D
      such that |D| <= Dmax.
      Non-fundamental positions need not contain the class number.
      Precisely, h [|D|] counts the number of quadratic forms [A, B, C]
      such that 0 <= |B| <= A <= C of discriminant D = B^2-4*A*C,
      and B>=0 if |B| = A or A = C. We may include in this count some or
      all of the non-primitive forms, since they belong to non-fundamental
      discriminants. */
{
   uint_cl_t Dmax2, A, B, D2, Amax, Alocmax, i;

   Dmax2 = Dmax / 2;
   for (i = 0; i <= Dmax2; i++)
      h [i] = 0;

   Amax = (uint_cl_t) sqrt (Dmax / 3.0);

   /* Consider forms with B=0. */
   for (A = 1; A <= Amax; A++)
      /* C = A, occurs as primitive form only for |D| = 4  */
      for (D2 = 2*A*A; D2 <= Dmax2; h [D2]++, D2 += 2*A);

   /* Consider forms with 0 < |B| = A <= C. */
   for (A = 1; A <= Amax; A++)
      /* C = A, occurs as primitive form only for |D| = 3 */
      for (D2 = 3 * A * A / 2; D2 <= Dmax2; h [D2]++, D2 += 2*A);

   /* Consider forms with 0 < |B| < A <= C. */
   for (B = 1; B < Amax; B++) {
      Alocmax = (uint_cl_t) sqrt ((Dmax + B * B) / 4);
      for (A = B + 1; A <= Alocmax; A++) {
         /* C = A */
         D2 = (2 * A - B) * (2 * A + B) / 2;
         h [D2]++;
         /* C = A + 1 */
         for (D2 += 2*A; D2 <= Dmax2; h [D2] += 2, D2 += 2*A);
      }
   }
}

/*****************************************************************************/

static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   int no_old, int no_new)
   /* Compute and return via qstar a list of "signed primes" suitable
      for dividing the discriminant to prove the primality of p, and via
      root their square roots modulo p. The entries in qstar are ordered
      by increasing absolute value.
      To make repeated calls possible for increasing the number of primes,
      no_old indicates the number of elements already present before the
      call, no_new the desired total number of elements.
      qstar and root need to be initialised to the correct size,
      and the entries of root need to be initialised as well. */
{
   long int q;
   unsigned int e;
   mpz_t r, z;

   mpz_init (r);
   mpz_init (z);
   e = cm_nt_mpz_tonelli_generator (r, z, p);

   q = (no_old == 0 ? 0 : qstar [no_old - 1]);

   while (no_old < no_new) {
      if (q == 0)
         q = -3;
      else if (q == -3)
         q = -4;
      else if (q == -4)
         q = 5;
      else if (q == 5)
         q = -7;
      else if (q == -7)
         q = -8;
      else if (q == -8)
         q = 8;
      else if (q == 8)
         q = -11;
      else {
         if (q > 0)
            q = cm_nt_next_prime (q);
         else
            q = cm_nt_next_prime (-q);
         if (q % 4 == 3)
            q = -q;
      }

      if (mpz_si_kronecker (q, p) == 1) {
         qstar [no_old] = q;
         cm_counter1++;
         cm_nt_mpz_tonelli_si_with_generator (root [no_old], q, p, e, r, z);
         no_old++;
      }
   }

   mpz_clear (r);
   mpz_clear (z);
}

/*****************************************************************************/

static int_cl_t** compute_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar, uint_cl_t Dmax, uint_cl_t hmaxprime,
   uint_cl_t *h)
   /* Given an array of no_qstar "signed primes" qstar (ordered by
      increasing absolute value), return an array of negative fundamental
      discriminants with factors from the list and of absolute value
      bounded above by Dmax, and return their number in no_d.
      Moreover, each discriminant must contain at least one prime larger
      in absolute value than the first no_qstar_old ones.
      If it is different from 0, then hmaxprime furthermore is an upper
      bound on the largest prime factor of the class number. Specifying it
      will result in the class numbers being factored, which takes
      additional time.
      Each element of the array is again an array (of fixed length)
      recording additional information on the discriminant. So far:
      0: discriminant
      1: class number h
      2: h/g, the class number relative to the genus field
      3: largest prime factor of h (only computed when hmax>0, otherwise
         it is set to 0)
      The additional input h is a precomputed array in which h [(-d)/2]
      contains the class number of the fundamental discriminant d.
      The task feels like it could be written in a few lines in a
      functional programming language; as a graph traversal, it is also
      readily implemented in C with "manual" backtracking. */
{
   int_cl_t **d;
      /* result */
   int_cl_t D;
      /* currently considered discriminant */
   int Dq [16];
      /* its prime factors as indices in qstar */
   int Dno;
      /* number of its prime factors, which is also the level in the tree */
   int Dqmax;
      /* largest index in qstar of one of its prime factors */
   unsigned int no;
      /* number of found discriminants so far */
   int no_factors;
   mpz_t tmp, bin;
   uint_cl_t hprime;
      /* largest prime factor of the class number, or 0 if not computed */
   int k;

   /* The algorithm assumes that at least the primes in qstar are of
      suitable size; otherwise, forget the largest ones. */
   while (   (qstar [no_qstar-1] > 0
              && (uint_cl_t) (qstar [no_qstar-1]) > Dmax)
          || (qstar [no_qstar-1] < 0
              && (uint_cl_t) (-qstar [no_qstar-1]) > Dmax))
      no_qstar--;

   /* Compute the maximal number of factors that can fit under Dmax. */
   k = 0;
   mpz_init_set_ui (tmp, 1);
   while (mpz_cmp_ui (tmp, Dmax) <= 0) {
      if (qstar [k] > 0)
         mpz_mul_ui (tmp, tmp, (unsigned long int) (qstar [k]));
      else
         mpz_mul_ui (tmp, tmp, (unsigned long int) (-qstar [k]));
      k++;
   }
   no_factors = k-1;

   /* Compute an upper bound on the possible number of discriminants. */
   mpz_init (bin);
   mpz_set_ui (tmp, 0);
   for (k = 1; k <= no_factors; k++) {
      mpz_bin_uiui (bin, no_qstar, k);
      mpz_add (tmp, tmp, bin);
   }
   if (mpz_cmp_ui (tmp, Dmax / 2) > 0)
      no = Dmax / 2;
   else
      no = mpz_get_ui (tmp);
   mpz_clear (bin);
   mpz_clear (tmp);

   d = (int_cl_t **) malloc (no * sizeof (int_cl_t *));

   no = 0;
   D = 1;
   Dno = 0;
   Dqmax = -1;
   /* Loop until we reach the last discriminant; in our depth first tree
      traversal, this is the one with only one prime factor, which is the
      largest one possible. */
   while (Dno != 1 || Dqmax != no_qstar - 1) {
      if (Dqmax < no_qstar - 1
          && (   (D > 0 && (uint_cl_t)   D  < Dmax)
              || (D < 0 && (uint_cl_t) (-D) < Dmax)))
         /* Add a level. */
         Dno++;
      else {
         /* Backtrack: If possible, stay at the same level and remove the
            current prime to replace it by the next candidate. If the
            current prime is already the largest one, we need to go one
            level up. Also if the current discriminant is already too
            large, there is no point in trying even larger primes, and we
            need to go one level up. */
         if (Dqmax == no_qstar - 1
             || (D > 0 && (uint_cl_t)   D  >= Dmax)
             || (D < 0 && (uint_cl_t) (-D) >= Dmax)) {
            D /= qstar [Dqmax];
            Dno--;
            Dqmax = Dq [Dno - 1];
         }
         /* On this level, remove the current prime. */
         D /= qstar [Dqmax];
      }
      /* Add one larger prime. After the "if" case above, Dqmax is the
         index of the currently largest prime; after the "else" case, it
         is the index of the prime we just removed. In both cases, it
         simply needs to be incremented. */
      Dqmax++;
      Dq [Dno - 1] = Dqmax;
      D *= qstar [Dqmax];

      /* Register the new discriminant if it satisfies the conditions. */
      if (D < 0
          && D % 16 != 0 /* only one of -4, -8 and 8 is included */
          && Dqmax >= no_qstar_old
          && (uint_cl_t) (-D) <= Dmax) {
         hprime = (hmaxprime > 0 ? cm_nt_largest_factor (h [(-D) / 2]) : 0);
         if (hprime <= hmaxprime) {
            d [no] = (int_cl_t *) malloc (4 * sizeof (int_cl_t));
            d [no][0] = D;
            d [no][1] = h [(-D) / 2];
            d [no][2] = d [no][1] >> (Dno - 1);
            d [no][3] = hprime;
            no++;
         }
      }
   }

   *no_d = no;
   d = realloc (d, no * sizeof (int_cl_t *));

   return d;
}

/*****************************************************************************/

static int disc_cmp (const void* d1, const void* d2)
{
   int_cl_t *D1, *D2;

   D1 = (*((int_cl_t **) d1));
   D2 = (*((int_cl_t **) d2));

   /* First sort by increasing h/g. */
   if (D1 [2] < D2 [2])
      return -1;
   else if (D1 [2] > D2 [2])
      return +1;
   else
      /* Then sort by increasing h. */
      if (D1 [1] < D2 [1])
         return -1;
      else if (D1 [1] > D2 [1])
         return +1;
      else
         /* Finally sort by increasing |d|. */
         if (D1 [0] < D2 [0])
            return +1;
         else if (D1 [0] > D2 [0])
            return -1;
         else
            return 0;
}

/*****************************************************************************/

static void trial_div (mpz_ptr l, mpz_srcptr n, mpz_srcptr primorialB)
   /* primorialB is supposed to be the product of the primes up to the
      smoothness bound B.
      The function removes all occurrences of the primes up to B from n
      and stores the result in l. */
{
   mpz_t mod, gcd;

   mpz_init (mod);
   mpz_init (gcd);

   mpz_set (l, n);
   mpz_mod (mod, primorialB, l);
   mpz_gcd (gcd, l, mod);
   while (mpz_cmp_ui (gcd, 1)) {
      mpz_divexact (l, l, gcd);
      mpz_gcd (gcd, l, mod);
   }

   mpz_clear (mod);
   mpz_clear (gcd);
}

/*****************************************************************************/

static bool is_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   mpz_srcptr root, const unsigned int delta, mpz_srcptr primorialB,
   int_cl_t d)
   /* The function tests whether d is a suitable discriminant to perform
      one step in the ECPP downrun from the (probable) prime N>=787 and
      returns the result. If true, n becomes the cardinality of the
      elliptic curve and l its largest prime factor; otherwise n and l are
      unchanged.
      root is a square root of d modulo N.
      delta >= 1 is the minimum number of bits to be gained in this step.
      primorialB is related to the trial division bound and is passed on
      to trial_div. */

{
   mpz_t t, v, co;
   mpz_ptr V;
   int twists;
   mpz_t card [6];
      /* The number of twists and the 2, 4, or 6 possible cardinalities
         of the curves. */
   size_t size_co, size_N;
   int i;
   bool ok;

   mpz_init (t);
   ok = false;
   if (d == -3 || d == -4) {
      mpz_init (v);
      V = v;
      if (d == -3)
         twists = 6;
      else
         twists = 4;
   }
   else {
      V = NULL;
      twists = 2;
   }
   if (cm_nt_mpz_cornacchia (t, V, N, root, d)) {

      /* Compute the cardinalities of all the twists. */
      mpz_init (co);
      for (i = 0; i < twists; i++)
         mpz_init (card [i]);
      mpz_add_ui (card [0], N, 1);
      mpz_add (card [1], card [0], t);
      if (d == -3) {
         /* The sextic twists have trace (\pm t \pm 3*v)/2. */
         mpz_mul_ui (v, v, 3);
         mpz_sub (v, v, t);
         mpz_divexact_ui (v, v, 2);
         mpz_add (card [2], card [0], v);
         mpz_sub (card [3], card [0], v);
         mpz_add (v, v, t);
         mpz_add (card [4], card [0], v);
         mpz_sub (card [5], card [0], v);
      }
      else if (d == -4) {
         /* The quartic twists have trace \pm 2*v. */
         mpz_mul_2exp (v, v, 1);
         mpz_sub (card [2], card [0], v);
         mpz_add (card [3], card [0], v);
      }
      mpz_sub (card [0], card [0], t);

      /* Cycle through the twists and look for a suitable cardinality. */
      for (i = 0; !ok && i < twists; i++) {
         cm_timer_continue (cm_timer5);
         cm_counter2++;
         trial_div (co, card [i], primorialB);
         cm_timer_stop (cm_timer5);
         /* We need to check whether co > (N^1/4 + 1)^2.
            Let N have e bits, that is, 2^(e-1) <= N < 2^e,
            and co have f bits, that is, 2^(f-1) <= co < 2^f. Then
            (N^(1/4) + 1)^2 = N^(1/2) * (1+1/N^(1/4))^2
                            < 2^(e/2) * sqrt (2) for N >= 781
            So it is sufficient to check that
            f-1 >= (e+1)/2, or equivalently f >= floor (e/2) + 2.
            We also want to gain at least one bit; otherwise we might
            even increase the prime a little bit. */
         cm_timer_continue (cm_timer3);
         size_co = mpz_sizeinbase (co, 2);
         size_N = mpz_sizeinbase (N, 2);
         if (size_co <= size_N - delta && size_co >= size_N / 2 + 2) {
            cm_counter3++;
            if (cm_nt_is_prime (co)) {
               ok = true;
               mpz_set (n, card [i]);
               mpz_set (l, co);
            }
         }
         cm_timer_stop (cm_timer3);
      }
      mpz_clear (co);
      for (i = 0; i < twists; i++)
         mpz_clear (card [i]);
   }
   if (d == -3 || d == -4)
      mpz_clear (v);
   mpz_clear (t);

   return ok;
}

/*****************************************************************************/

static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const unsigned int delta, mpz_srcptr primorialB)
   /* Given a (probable) prime N>=787, return a suitable CM discriminant
      and return the cardinality of an associated elliptic curve in n and
      its largest prime factor in l.
      Dmax, hmaxprime and h are passed through to compute_discriminants.
      delta >= 1 is passed through as the minimum number of bits to be
      gained in this step.
      primorialB is passed through to trial division. */
{
   int no_qstar_old, no_qstar;
   long int qstar [1000];
   mpz_t root [1000];
   int i, j;
   int_cl_t d;
   int_cl_t **dlist;
   int no_d;
   mpz_t Droot;
   bool ok;

   ok = false;
   no_qstar_old = 0;
   no_qstar = 10;
   mpz_init (Droot);

   do {
      /* Extend the prime and square root list. */
      cm_timer_continue (cm_timer1);
      for (i = no_qstar_old; i < no_qstar; i++)
         mpz_init (root [i]);
      compute_qstar (qstar, root, N, no_qstar_old, no_qstar);
      cm_timer_stop (cm_timer1);

      /* Precompute a list of potential discriminants for fastECPP. */
      cm_timer_continue (cm_timer2);
      dlist = compute_discriminants (&no_d, qstar, no_qstar_old, no_qstar,
         Dmax, hmaxprime, h);
      qsort (dlist, no_d, sizeof (int_cl_t *), disc_cmp);
      cm_timer_stop (cm_timer2);

      /* Search for the first suitable discriminant in the list. */
      i = -1;
      do {
         i++;
         d = dlist [i][0];
         /* Compute the square root of the discriminant by multiplying the
            roots of all its prime factors. d can be trial divided over qstar,
            but the "even primes" need special care, in particular for the
            sign of 8. */
         mpz_set_ui (Droot, 1);
         for (j = 0; j < no_qstar; j++)
            if (qstar [j] % 2 != 0 && d % qstar [j] == 0) {
               mpz_mul (Droot, Droot, root [j]);
               mpz_mod (Droot, Droot, N);
               d /= qstar [j];
            }
         if (d != 1)
            for (j = 0; j < no_qstar; j++)
               if (d == qstar [j]) {
                  mpz_mul (Droot, Droot, root [j]);
                  mpz_mod (Droot, Droot, N);
               }
         d = dlist [i][0];
         ok = is_ecpp_discriminant (n, l, N, Droot, delta, primorialB, d);
      } while (!ok && i < no_d - 1);

      if (!ok) {
         no_qstar_old = no_qstar;
         if (no_qstar <= 10)
            no_qstar += 10;
         else
            no_qstar += 20;
         if (no_qstar > 1000) {
            printf ("Error in find_ecpp_discriminant: qstar array too short\n");
            exit (1);
         }
      }

      /* Free the discriminant list. */
      for (i = 0; i < no_d; i++)
         free (dlist [i]);
      free (dlist);

   } while (!ok);

   mpz_clear (Droot);
   for (i = 0; i < no_qstar; i++)
      mpz_clear (root [i]);

   return d;
}

/*****************************************************************************/

mpz_t** cm_ecpp1 (int *depth, mpz_srcptr p, bool verbose, bool debug)
   /* Compute the first step of the ECPP certificate; this is the downrun
      part with the parameters of the elliptic curves.
      The return value is a newly allocated array of depth entries, each
      of which is an array of length 4, containing in this order
      - p_i, a prime to be certified;
      - d_i, the discriminant;
      - n_i, the cardinality of the elliptic curve;
      - l_i, the prime order dividing this cardinality.
      The downrun stops as soon as the prime is less than 2^64. */

{
   const unsigned long int B = 1000000;
   const unsigned int delta = (unsigned int) (log2 (B) / 2) + 1;
      /* According to [FrKlMoWi04] the average factor removed by trial
         division up to B, assuming that what remains is prime, is B;
         we impose half of this number of bits as the minimal gain. */
   const uint_cl_t hmaxprime = 30;
   mpz_t N, primorialB;
   mpz_t** c;
   uint_cl_t *h;
   uint_cl_t Dmax;
   int_cl_t d;
   cm_timer clock;
   int i;

   /* Precompute class numbers. */
   cm_timer_start (clock);
   Dmax = mpz_sizeinbase (p, 2);
   Dmax = ((Dmax * Dmax) >> 3) << 1;
   h = (uint_cl_t *) malloc ((Dmax / 2 + 1) * sizeof (uint_cl_t));
   compute_h (h, Dmax);
   cm_timer_stop (clock);
   if (verbose)
      printf ("-- Time for class numbers up to Dmax=%"PRIucl": %5.1f\n",
         Dmax, cm_timer_get (clock));

   /* Precompute primorial for trial division. */
   cm_timer_start (clock);
   mpz_init (primorialB);
   mpz_primorial_ui (primorialB, B);
   cm_timer_stop (clock);
   if (verbose)
      printf ("-- Time for primorial: %5.1f\n", cm_timer_get (clock));

   mpz_init_set (N, p);
   *depth = 0;
   c = (mpz_t**) malloc (*depth);
   cm_timer_reset (cm_timer1);
   cm_timer_reset (cm_timer2);
   cm_timer_reset (cm_timer3);
   cm_timer_reset (cm_timer4);
   cm_timer_reset (cm_timer5);
   while (mpz_sizeinbase (N, 2) > 64) {
      c = (mpz_t**) realloc (c, (*depth + 1) * sizeof (mpz_t *));
      c [*depth] = (mpz_t *) malloc (4 * sizeof (mpz_t));
      for (i = 0; i < 4; i++)
         mpz_init (c [*depth][i]);
      mpz_set (c [*depth][0], N);
      cm_counter1 = 0;
      cm_counter2 = 0;
      cm_counter3 = 0;
      cm_timer_start (clock);
      d = find_ecpp_discriminant (c [*depth][2], c [*depth][3], N, Dmax,
         hmaxprime, h, delta, primorialB);
      cm_timer_stop (clock);
      if (verbose) {
         printf ("-- Time for discriminant %8"PRIicl" for %4li bits: %5.1f\n",
            d, mpz_sizeinbase (N, 2), cm_timer_get (clock));
         if (debug) {
            printf ("%5i qstar: %.1f, disclist: %.1f\n",
                  cm_counter1, cm_timer_get (cm_timer1), cm_timer_get (cm_timer2));
            printf ("%5i Trial div: %.1f\n", cm_counter2,
                  cm_timer_get (cm_timer5));
            printf ("%5i is_prime: %.1f\n", cm_counter3,
                  cm_timer_get (cm_timer3));
         }
      }
      mpz_set_si (c [*depth][1], d);
      mpz_set (N, c [*depth][3]);
      (*depth)++;
   }

   free (h);
   mpz_clear (primorialB);
   mpz_clear (N);

   return c;
}

/*****************************************************************************/

void cm_ecpp (mpz_srcptr N, const char* modpoldir, bool pari, bool tower,
   bool print, bool verbose, bool debug)
   /* Assuming that N is a (probable) prime, compute an ECPP certificate.
      modpoldir gives the directory where modular polynomials are stored;
      it is passed through to the function computing a curve from a root
      of the class polynomial.
      pari indicates whether the first, downrun step of PARI is used
      instead of the built-in function.
      tower indicates whether a class field tower decomposition is used
      instead of only the class polynomial.
      print indicates whether the result is printed.
      verbose indicates whether intermediate computations output progress
      information.
      debug indicates whether additional developer information (mainly
      timings and counters for tuning) is output; this is done only in
      the case that verbose is set as well. */

{
   int depth;
   mpz_t **cert;
   mpz_srcptr p, n, l;
   int_cl_t d;
   mpz_t t, co, a, b, x, y;
   cm_param_t param, new_param;
   double hf, new_hf;
   cm_class_t c;
   int i, j;
   cm_timer clock, clock2, clock3, clock4, clock5;

   mpz_init (t);
   mpz_init (co);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);

   cm_timer_start (clock);
   cm_timer_start (clock2);
   if (pari)
      cert = cm_pari_ecpp1 (&depth, N);
   else
      cert = cm_ecpp1 (&depth, N, verbose, debug);
   cm_timer_stop (clock2);
   if (verbose)
      printf ("--- Time for first ECPP step, depth %i:  %.1f\n", depth,
         cm_timer_get (clock2));

   cm_timer_start (clock2);
   cm_timer_reset (clock4);
   cm_timer_reset (clock5);
   cm_timer_reset (cm_timer1);
   cm_timer_reset (cm_timer2);
   cm_timer_reset (cm_timer3);
   cm_timer_reset (cm_timer4);
   cm_timer_reset (cm_timer5);
   cm_timer_reset (cm_timer6);
   if (print)
      printf ("c = [");
   for (i = 0; i < depth; i++) {
      cm_timer_start (clock3);
      p = cert [i][0];
      d = mpz_get_si (cert [i][1]);
      n = cert [i][2];
      l = cert [i][3];

      mpz_add_ui (t, p, 1);
      mpz_sub (t, t, n);
      mpz_divexact (co, n, l);

      cm_timer_continue (clock4);
      /* Find the invariant with the largest height factor from those
         not requiring modular polynomials for retrieving the curve. */
      /* First test the non-parametric invariants in the good order. */
      if (   !cm_param_init (param, d, CM_INVARIANT_WEBER, false)
          && !(((d - 1) % 8 == 0
               && cm_param_init (param, 4*d, CM_INVARIANT_WEBER, false)))
          && !cm_param_init (param, d, CM_INVARIANT_GAMMA2, false)
          && !cm_param_init (param, d, CM_INVARIANT_GAMMA3, false))
          cm_param_init (param, d, CM_INVARIANT_J, true);
      hf = cm_class_height_factor (param);
      /* Atkin invariants have excellent factors between 24 and 36, but
         the Hecke operators are slow to compute. So do not use them.
         Simple eta is currently hard-wired to w_3^12 with a factor of 4,
         which is still better than all others currently available when
         Weber does not work. */
      if (cm_param_init (new_param, d, CM_INVARIANT_SIMPLEETA, false)) {
         new_hf = cm_class_height_factor (new_param);
         if (new_hf > hf) {
            param [0] = new_param [0];
            hf = new_hf;
         }
      }
      /* FIXME: For double and multiple eta quotients, one needs to take
         the degree of the modular polynomial into account. */

      /* Compute one of the class field tower or the class polynomial. */
      cm_class_init (c, param, false);
      cm_class_compute (c, param, !tower, tower, false);
      cm_timer_stop (clock4);

      cm_timer_continue (clock5);
      cm_curve_and_point (a, b, x, y, param, c, p, l, co,
         modpoldir, false);
      cm_class_clear (c);
      cm_timer_stop (clock5);
      cm_timer_stop (clock3);

      if (verbose) {
         printf ("-- Time for discriminant %8"PRIicl" with invariant %c "
            "for %4li bits: %5.1f\n",
            d, param->invariant, mpz_sizeinbase (p, 2),
            cm_timer_get (clock3));
         if (debug) {
            printf ("   CM:    %5.1f\n", cm_timer_get (clock4));
            printf ("   roots: %5.1f\n", cm_timer_get (cm_timer1));
            printf ("   curve: %5.1f\n", cm_timer_get (cm_timer2));
            printf ("     random:   %5.1f\n", cm_timer_get (cm_timer3));
            printf ("     multiply: %5.1f\n", cm_timer_get (cm_timer4));
            printf ("       dbl: %5.1f\n", cm_timer_get (cm_timer5));
            printf ("       add: %5.1f\n", cm_timer_get (cm_timer6));
         }
      }

      if (print) {
         printf ("[");
         mpz_out_str (stdout, 10, p);
         printf (", ");
         mpz_out_str (stdout, 10, t);
         printf (", ");
         mpz_out_str (stdout, 10, co);
         printf (", ");
         mpz_out_str (stdout, 10, a);
         printf (", [");
         mpz_out_str (stdout, 10, x);
         printf (", ");
         mpz_out_str (stdout, 10, y);
         printf ("]]");
         if (i != depth - 1)
            printf (", ");
      }
   }
   if (print)
      printf ("];\n");

   mpz_clear (t);
   mpz_clear (co);
   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (x);
   mpz_clear (y);
   for (i = 0; i < depth; i++) {
      for (j = 0; j < 4; j++)
         mpz_clear (cert [i][j]);
      free (cert [i]);
   }
   free (cert);
   cm_timer_stop (clock2);
   cm_timer_stop (clock);
   if (verbose) {
      printf ("--- Time for CM:               %.1f\n", cm_timer_get (clock4));
      printf ("--- Time for curves:           %.1f\n", cm_timer_get (clock5));
      printf ("--- Time for second ECPP step: %.1f\n", cm_timer_get (clock2));
      printf ("--- Total time for ECPP:       %.1f\n", cm_timer_get (clock));
   }
}

/*****************************************************************************/
/*****************************************************************************/
