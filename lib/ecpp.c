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

static bool is_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   mpz_srcptr root, int_cl_t d);
static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N);
static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   int no);
static int_cl_t* compute_discriminants (int *no_d, long int *qstar,
   int no_qstar, int no_factors, uint_cl_t Dmax);
static int disc_cmp (const void* d1, const void* d2);

/*****************************************************************************/

static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   int no)
   /* Compute and return via qstar a list of "signed primes" suitable
      for dividing the discriminant to prove the primality of p, and via
      root their square roots modulo p. The entries in qstar are ordered
      by increasing absolute value.
      no indicates the desired number of elements.
      qstar and root need to be initialised to the correct size,
      and the entries of root need to be initialised as well. */
{
   int i;
   long int q;

   i = 0;
   q = -3;
   while (i < no) {
      if (mpz_si_kronecker (q, p) == 1) {
         qstar [i] = q;
         cm_nt_mpz_tonelli (root [i], q, p);
         i++;
      }
      if (q == -3)
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
   }
}

/*****************************************************************************/

static int_cl_t* compute_discriminants (int *no_d, long int *qstar,
   int no_qstar, int no_factors, uint_cl_t Dmax)
   /* Given an array of no_qstar "signed primes" qstar (ordered by
      increasing absolute value), return an array of negative fundamental
      discriminants with at most no_factors prime factors from the list,
      and return their number in no_d.
      If it is different from 0, then Dmax furthermore is an upper bound
      on the absolute value of the discriminant.
      The task feels like it could be written in a few lines in a
      functional programming language; as a graph traversal, it is also
      readily implemented in C with "manual" backtracking. */
{
   int_cl_t *d;
      /* result */
   int_cl_t D;
      /* currently considered discriminant */
   int Dq [16];
      /* its prime factors as indices in qstar */
   int Dno;
      /* number of its prime factors, which is also the level in the tree */
   int Dqmax;
      /* largest index in qstar of one of its prime factors */
   int no;
      /* number of found discriminants so far */
   mpz_t Dz;
   int k;

   if (Dmax != 0)
      /* The algorithm assumes that at least the primes in qstar are of
         suitable size; otherwise, forget the largest ones. */
      while (   (qstar [no_qstar-1] > 0
                 && (uint_cl_t) (qstar [no_qstar-1]) > Dmax)
             || (qstar [no_qstar-1] < 0
                 && (uint_cl_t) (-qstar [no_qstar-1]) > Dmax))
         no_qstar--;

   if (no_factors > no_qstar)
      no_factors = no_qstar;

   if (Dmax != 0) {
      /* no_factors may be larger than the actual number of factors
         that can fit under Dmax; since the elements of qstar are
         sorted, the product of the first no_factors elements must
         fit. */
      k = 0;
      mpz_init_set_ui (Dz, 1);
      while (k < no_factors && mpz_cmp_ui (Dz, Dmax) <= 0) {
         if (qstar [k] > 0)
            mpz_mul_ui (Dz, Dz, (unsigned long int) (qstar [k]));
         else
            mpz_mul_ui (Dz, Dz, (unsigned long int) (-qstar [k]));
         k++;
      }
      if (mpz_cmp_ui (Dz, Dmax) > 0)
         no_factors = k-1;
      mpz_clear (Dz);
   }

   /* Compute an upper bound on the possible number of discriminants. */
   no = 0;
   for (k = 1; k <= no_factors; k++)
      no += cm_nt_binomial (no_qstar, k);
   d = (int_cl_t *) malloc (no * sizeof (int_cl_t));

   no = 0;
   D = 1;
   Dno = 0;
   Dqmax = -1;

   /* Loop until we reach the last discriminant; in our depth first tree
      traversal, this is the one with only one prime factor, which is the
      largest one possible. */
   while (Dno != 1 || Dqmax != no_qstar - 1) {
      if (Dno < no_factors && Dqmax < no_qstar - 1)
         /* Add a level. */
         Dno++;
      else {
         /* Backtrack: Find the largest previous level where Dqmax
            can be increased. */
         while (Dqmax == no_qstar - 1) {
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
      if (D < 0 && (Dmax == 0 || (uint_cl_t) (-D) <= Dmax)) {
         d [no] = D;
         no++;
      }
   }

   *no_d = no;
   d = realloc (d, no * sizeof (int_cl_t));

   return d;
}

/*****************************************************************************/

static int disc_cmp (const void* d1, const void* d2)
{
   int_cl_t D1, D2;

   D1 = *((int_cl_t *) d1);
   D2 = *((int_cl_t *) d2);

   return (D1 < D2 ? +1 : (D1 == D2 ? 0 : -1));
}

/*****************************************************************************/

static bool is_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   mpz_srcptr root, int_cl_t d)
   /* The function tests whether d is a suitable discriminant to perform
      one step in the ECPP downrun from the (probable) prime N>=787 and
      returns the result. If true, n becomes the cardinality of the
      elliptic curve and l its largest prime factor; otherwise n and l are
      unchanged.
      root is a square root of d modulo N. */

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
         cm_pari_trialdiv (co, card [i], 1000000);
         /* We need to check whether co > (N^1/4 + 1)^2.
            Let N have e bits, that is, 2^(e-1) <= N < 2^e,
            and co have f bits, that is, 2^(f-1) <= co < 2^f. Then
            (N^(1/4) + 1)^2 = N^(1/2) * (1+1/N^(1/4))^2
                            < 2^(e/2) * sqrt (2) for N >= 781
            So it is sufficient to check that
            f-1 >= (e+1)/2, or equivalently f >= floor (e/2) + 2.
            We also want to gain at least one bit; otherwise we might
            even increase the prime a little bit. */
         size_co = mpz_sizeinbase (co, 2);
         size_N = mpz_sizeinbase (N, 2);
         if (size_co < size_N && size_co >= size_N / 2 + 2
            && cm_nt_is_prime (co)) {
            ok = true;
            mpz_set (n, card [i]);
            mpz_set (l, co);
         }
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

static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N)
   /* Given a (probable) prime N>=787, return a suitable CM discriminant
      and return the cardinality of an associated elliptic curve in n and
      its largest prime factor in l. */
{
   int no_qstar = 20;
   int no_factors = 3;
   uint_cl_t Dmax = 0;
   long int *qstar;
   mpz_t *root;
   int i, j;
   int_cl_t d;
   int_cl_t *dlist;
   int no_d;
   mpz_t Droot;
   bool ok;

   /* Prepare the prime and square root list. */
   qstar = (long int *) malloc (no_qstar * sizeof (long int));
   root = (mpz_t *) malloc (no_qstar * sizeof (mpz_t));
   for (i = 0; i < no_qstar; i++)
      mpz_init (root [i]);
   compute_qstar (qstar, root, N, no_qstar);

   /* Precompute a list of potential discriminants for fastECPP. */
   dlist = compute_discriminants (&no_d, qstar, no_qstar, no_factors, Dmax);
   qsort (dlist, no_d, sizeof (int_cl_t), disc_cmp);

   /* Search for the first suitable discriminant in the list. */
   mpz_init (Droot);
   i = -1;
   ok = false;
   do {
      i++;
      d = dlist [i];
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
      d = dlist [i];
      ok = is_ecpp_discriminant (n, l, N, Droot, d);
   } while (!ok && i < no_d - 1);
   if (!ok) {
      printf ("***** Error in find_ecpp_discriminant: qstar too short\n");
      exit (1);
   }
   mpz_clear (Droot);

   free (dlist);
   for (i = 0; i < no_qstar; i++)
      mpz_clear (root [i]);
   free (root);
   free (qstar);

   return d;
}


/*****************************************************************************/

mpz_t** cm_ecpp1 (int *depth, mpz_srcptr p, bool verbose)
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
   mpz_t N;
   mpz_t** c;
   int_cl_t d;
   cm_timer clock;
   int i;

   mpz_init_set (N, p);
   *depth = 0;
   c = (mpz_t**) malloc (*depth);
   while (mpz_sizeinbase (N, 2) > 64) {
      c = (mpz_t**) realloc (c, (*depth + 1) * sizeof (mpz_t *));
      c [*depth] = (mpz_t *) malloc (4 * sizeof (mpz_t));
      for (i = 0; i < 4; i++)
         mpz_init (c [*depth][i]);
      mpz_set (c [*depth][0], N);
      cm_timer_start (clock);
      d = find_ecpp_discriminant (c [*depth][2], c [*depth][3], N);
      cm_timer_stop (clock);
      if (verbose)
         printf ("-- Time for discriminant %6"PRIicl" for %4li bits: %5.1f\n",
            d, mpz_sizeinbase (N, 2), cm_timer_get (clock));
      mpz_set_si (c [*depth][1], d);
      mpz_set (N, c [*depth][3]);
      (*depth)++;
   }

   mpz_clear (N);

   return c;
}

/*****************************************************************************/

void cm_ecpp (mpz_srcptr N, const char* modpoldir, bool pari, bool tower,
   bool print, bool verbose)
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
      information. */
{
   int depth;
   mpz_t **cert;
   mpz_srcptr p, n, l;
   int_cl_t d;
   mpz_t t, co, a, b, x, y;
   cm_param_t param;
   cm_class_t c;
   int i, j;
   cm_timer clock, clock1, clock2;

   mpz_init (t);
   mpz_init (co);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);

   cm_timer_start (clock);
   if (pari)
      cert = cm_pari_ecpp1 (&depth, N);
   else
      cert = cm_ecpp1 (&depth, N, verbose);
   cm_timer_stop (clock);
   if (verbose)
      printf ("--- Time for first ECPP step:  %.1f\n", cm_timer_get (clock));

   cm_timer_start (clock);
   cm_timer_reset (clock1);
   cm_timer_reset (clock2);
   if (print)
      printf ("c = [");
   for (i = 0; i < depth; i++) {
      p = cert [i][0];
      d = mpz_get_si (cert [i][1]);
      n = cert [i][2];
      l = cert [i][3];

      mpz_add_ui (t, p, 1);
      mpz_sub (t, t, n);
      mpz_divexact (co, n, l);

      cm_timer_continue (clock1);
      if (d % 3 != 0)
         cm_param_init (param, d, CM_INVARIANT_GAMMA2, false);
      else if (d % 2 != 0)
         cm_param_init (param, d, CM_INVARIANT_GAMMA3, false);
      else
         cm_param_init (param, d, CM_INVARIANT_J, false);
      cm_class_init (c, param, false);
      cm_class_compute (c, param, !tower, tower, false);
      cm_timer_stop (clock1);
      cm_timer_continue (clock2);
      cm_curve_and_point (a, b, x, y, param, c, p, l, co,
         modpoldir, false);
      cm_timer_stop (clock2);
      cm_class_clear (c);

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
   cm_timer_stop (clock);
   if (verbose) {
      printf ("--- Time for CM:               %.1f\n", cm_timer_get (clock1));
      printf ("--- Time for curves:           %.1f\n", cm_timer_get (clock2));
      printf ("--- Time for second ECPP step: %.1f\n", cm_timer_get (clock));
   }
}

/*****************************************************************************/
/*****************************************************************************/
