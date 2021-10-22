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

static void compute_h_chunk (uint_cl_t *h, uint_cl_t Dmin, uint_cl_t Dmax);
static void compute_h (uint_cl_t *h, uint_cl_t Dmax);
static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   long int *q, int no);
static int_cl_t** compute_signed_discriminants (int *no_d, long int *qstar,
   int no_qstar, uint_cl_t Dmax, int sign);
static int_cl_t** compute_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h);
static int_cl_t* compute_sorted_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h);
static int disc_cmp (const void* d1, const void* d2);
static void mpz_tree_mod (mpz_t *mod, mpz_srcptr n, mpz_t *m, int no_m);
static void trial_div_batch (mpz_t *l, mpz_t *n, int no_n,
   mpz_srcptr primorialB);
static int curve_cardinalities (mpz_t *n, mpz_srcptr N, mpz_srcptr root,
   int_cl_t d);
static mpz_t* compute_cardinalities (int *no_card, int_cl_t **card_d,
   mpz_srcptr N, int no_d, mpz_t *root, int_cl_t *d);
static int_cl_t contains_ecpp_discriminant (mpz_ptr n, mpz_ptr l,
   mpz_srcptr N, const unsigned int delta, mpz_srcptr primorialB,
   int_cl_t *d, int no_d, mpz_t *root);
static void root_of_d (mpz_t *Droot, int_cl_t *d, int no_d, mpz_srcptr N,
   long int *qstar, int no_qstar, mpz_t *root);
static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const unsigned int delta, mpz_srcptr primorialB,
   bool verbose, bool debug);
static void ecpp_param_init (cm_param_ptr param, uint_cl_t d);
static void ecpp2_one_step (mpz_t *cert2, mpz_t *cert1,
   const char* modpoldir, bool tower, bool verbose, bool debug);
static void cm_ecpp2 (mpz_t **cert2, mpz_t **cert1, int depth,
   const char* modpoldir, bool tower, bool verbose, bool debug);

/*****************************************************************************/

static void compute_h_chunk (uint_cl_t *h, uint_cl_t Dmin, uint_cl_t Dmax)
   /* Assuming that h is an array of length (Dmax-Dmin)/2 for Dmin and Dmax
      both divisible by 4, compute in h [(|D|-Dmin)/2-1] the class number
      for fundamental discriminants D such that Dmin < |D| <= Dmax.
      Non-fundamental positions need not contain the class number.
      Precisely, h [(|D|-Dmin)/2-1] counts the number of quadratic forms
      [A, B, C] such that 0 <= |B| <= A <= C of discriminant D = B^2-4*A*C,
      and B>=0 if |B| = A or A = C. We may include in this count some or
      all of the non-primitive forms, since they belong to non-fundamental
      discriminants. */
{
   uint_cl_t Dmin2, Dmax2, length, D2, A, B, A2, Amax, Alocmax, i;

   if (Dmin % 4 != 0 || Dmax % 4 != 0) {
      printf ("***** Error: compute_h_chunk called with parameters not "
              "divisible by 4.\n");
      exit (1);
   }

   Dmin2 = Dmin / 2;
   Dmax2 = Dmax / 2;
   length = Dmax2 - Dmin2;
   for (i = 0; i < length; i++)
      h [i] = 0;

   Amax = (uint_cl_t) sqrt (Dmax / 3.0);

   /* Consider forms with B=0. */
   Alocmax = (uint_cl_t) sqrt (Dmax2);
   for (A = 1; A <= Alocmax; A++) {
      A2 = 2 * A;
      /* Compute D2 = |D| / 2 corresponding to C = A. */
      D2 = A2 * A;
      if (D2 <= Dmin2)
         D2 += A2 * ((Dmin2 - D2 + A2) / A2);
      for (i = D2 - Dmin2 - 1; i < length; h [i]++, i += A2);
   }

   /* Consider forms with 0 < B = A <= C. */
   for (A = 1; A <= Amax; A++) {
      A2 = 2 * A;
      /* Compute D2 corresponding to C = A. */
      D2 = 3*A*A / 2;
      if (D2 <= Dmin2)
         D2 += A2 * ((Dmin2 - D2 + A2) / A2);
      for (i = D2 - Dmin2 - 1; i < length; h [i]++, i += A2);
   }

   /* Consider forms with 0 < |B| < A <= C. */
   for (B = 1; B < Amax; B++) {
      Alocmax = (uint_cl_t) sqrt ((Dmax + B * B) / 4.0);
      for (A = B + 1; A <= Alocmax; A++) {
         A2 = 2 * A;
         /* Compute D2 corresponding to C = A; if this is in the correct
            range, then the ambiguous form needs to be counted once. */
         D2 = (A2 - B) * (A2 + B) / 2;
         if (D2 > Dmin2) {
            h [D2 - Dmin2 - 1]++;
            D2 += A2;
         }
         else
            D2 += A2 * ((Dmin2 - D2 + A2) / A2);
         for (i = D2 - Dmin2 - 1; i < length; h [i] += 2, i += A2);
      }
   }
}

/*****************************************************************************/

static void compute_h (uint_cl_t *h, uint_cl_t Dmax)
   /* The function behaves as compute_h_chunk (h, 0, Dmax), but computing
      the class numbers in ranges of 100000 (that is, 50000 discriminants
      at a time). Experimentally, this optimises the running time on my
      laptop for Dmax=10^7. The optimal value probably depends on the cache
      size. */
{
   const uint_cl_t size = 100000;
   int chunks, i;

   if (Dmax % 4 != 0) {
      printf ("***** Error: compute_h called with parameter not "
              "divisible by 4.\n");
      exit (1);
   }

   chunks = (Dmax + size - 1) / size;
   for (i = 0; i < chunks - 1; i++)
      compute_h_chunk (h + i * size / 2, i * size, (i+1) * size);
   compute_h_chunk (h + i * size / 2, i * size, Dmax);
}

/*****************************************************************************/

static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   long int *q, int no)
   /* Compute and return via qstar an array of no "signed primes" suitable
      for dividing the discriminant to prove the primality of p, and via
      root their square roots modulo p. The entries in qstar are ordered
      by increasing absolute value. On first call, *q should be 0; it is
      then replaced by the largest found signed prime, and should be given
      for later calls to continue with the next signed primes.
      qstar and root need to be initialised to the correct size,
      and the entries of root need to be initialised as well. */
{
   unsigned int e;
   mpz_t r, z;
   int i;

   mpz_init (r);
   mpz_init (z);

   i = 0;
   while (i < no) {
      if (*q == 0)
         *q = -3;
      else if (*q == -3)
         *q = -4;
      else if (*q == -4)
         *q = 5;
      else if (*q == 5)
         *q = -7;
      else if (*q == -7)
         *q = -8;
      else if (*q == -8)
         *q = 8;
      else if (*q == 8)
         *q = -11;
      else {
         if (*q > 0)
            *q = cm_nt_next_prime (*q);
         else
            *q = cm_nt_next_prime (-*q);
         if (*q % 4 == 3)
            *q = -*q;
      }

      if (mpz_si_kronecker (*q, p) == 1) {
         qstar [i] = *q;
         i++;
      }
   }

   e = cm_nt_mpz_tonelli_generator (r, z, p);
   cm_counter1 += no;
   for (i = 0; i < no; i++)
      cm_nt_mpz_tonelli_si_with_generator (root [i], qstar [i], p, e, r, z);

   mpz_clear (r);
   mpz_clear (z);
}

/*****************************************************************************/

static int_cl_t** compute_signed_discriminants (int *no_d, long int *qstar,
   int no_qstar, uint_cl_t Dmax, int sign)
   /* Given an array of no_qstar "signed primes" qstar (ordered by
      increasing absolute value), return an array of fundamental
      discriminants of the given sign with factors from the list and of
      absolute value bounded above by Dmax, and return their number in no_d.
      For our purposes, 1 counts as a positive discriminant.
      Each element of the array is again an array (of fixed length)
      recording additional information on the discriminant. So far:
      0: discriminant
      1: number of prime factors (this is needed later for h/g)
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
   no ++;

   d = (int_cl_t **) malloc (no * sizeof (int_cl_t *));

   no = 0;
   D = 1;
   Dno = 0;
   Dqmax = -1;
   if (sign == 1) {
      d [no] = (int_cl_t *) malloc (2 * sizeof (int_cl_t));
      d [no][0] = 1;
      d [no][1] = 0;
      no++;
   }
   if (no_qstar >= 1) {
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
         if (D % 16 != 0 /* only one of -4, -8 and 8 is included */
             && (   (sign > 0 && D > 0 && (uint_cl_t)   D  <= Dmax)
                 || (sign < 0 && D < 0 && (uint_cl_t) (-D) <= Dmax))) {
            d [no] = (int_cl_t *) malloc (2 * sizeof (int_cl_t));
            d [no][0] = D;
            d [no][1] = Dno;
            no++;
         }
      }
   }

   *no_d = no;
   d = realloc (d, no * sizeof (int_cl_t *));

   return d;
}

/*****************************************************************************/

static int_cl_t** compute_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h)
   /* Given an array of no_qstar_old + no_qstar_new "signed primes" qstar
      (ordered by increasing absolute value), return an array of negative
      fundamental discriminants with factors from the list and of absolute
      value bounded above by Dmax, and return their number in no_d.
      Moreover, each discriminant must contain at least one prime larger
      from the last no_qstar_new ones.
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
      contains the class number of the fundamental discriminant d. */
{
   int_cl_t **d;
      /* result */
   int_cl_t ***d_part;
      /* partial results, one for each new element of qstar */
   int *no_part;
      /* lengths of the partial results */
   int no;
   int_cl_t qnew, D, Dno;
   uint_cl_t hprime;
      /* largest prime factor of the class number, or 0 if not computed */
   int i, j;

   d_part = (int_cl_t ***) malloc (no_qstar_new * sizeof (int_cl_t **));
   no_part = (int *) malloc (no_qstar_new * sizeof (int));
   for (i = 0; i < no_qstar_new; i++) {
      qnew = qstar [no_qstar_old + i];
      /* Compute all discriminants with primes less than qnew in absolute
         value that can be multiplied by qnew to obtain a suitable
         candidate. */
      if (qnew > 0)
         d_part [i] = compute_signed_discriminants (&(no_part [i]),
            qstar, no_qstar_old + i, Dmax / qnew, -1);
      else
         d_part [i] = compute_signed_discriminants (&(no_part [i]),
            qstar, no_qstar_old + i, Dmax / (-qnew), +1);
   }

   /* Filter all suitable discriminants and concatenate the lists. */
   no = 0;
   for (i = 0; i < no_qstar_new; i++)
      no += no_part [i];
   d = (int_cl_t **) malloc (no * sizeof (int_cl_t *));

   no = 0;
   for (i = 0; i < no_qstar_new; i++) {
      qnew = qstar [no_qstar_old + i];
      for (j = 0; j < no_part [i]; j++) {
         D = qnew * d_part [i][j][0];
         if (D % 16 != 0) {
            Dno = 1 + d_part [i][j][1];
            if (Dno <= max_factors) {
               hprime = (hmaxprime > 0 ?
                         cm_nt_largest_factor (h [(-D) / 2 - 1]) : 0);
               if (hprime <= hmaxprime) {
                  d [no] = (int_cl_t *) malloc (4 * sizeof (int_cl_t));
                  d [no][0] = D;
                  d [no][1] = h [(-D) / 2 - 1];
                  d [no][2] = d [no][1] >> (Dno - 1);
                  d [no][3] = hprime;
                  no++;
               }
            }
         }
         free (d_part [i][j]);
      }
      free (d_part [i]);
   }
   free (d_part);
   free (no_part);

   *no_d = no;
   d = realloc (d, no * sizeof (int_cl_t *));

   return d;
}

/*****************************************************************************/

static int_cl_t* compute_sorted_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h)
   /* The function takes the same parameters as compute_discriminants (and
      most of them are just passed through). But instead of an unsorted
      double array, it returns a simple array of discriminants sorted
      according to the usefulness in ECPP. */
{
   int_cl_t **dlist;
   int_cl_t *d;
   int i;

   dlist = compute_discriminants (no_d, qstar, no_qstar_old, no_qstar_new,
      max_factors, Dmax, hmaxprime, h);
   qsort (dlist, *no_d, sizeof (int_cl_t *), disc_cmp);

   d = (int_cl_t *) malloc (*no_d * sizeof (int_cl_t));
   for (i = 0; i < *no_d; i++) {
      d [i] = dlist [i][0];
      free (dlist [i]);
   }
   free (dlist);

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

static void mpz_tree_mod (mpz_t *mod, mpz_srcptr n, mpz_t *m, int no_m)
   /* Given a positive integer n and an array of no_m positive moduli in m,
      compute n mod all the m and return them in mod, which needs to
      provide sufficient space and initialised entries. For the
      implementation to be efficient, it is required that the product of
      all the m be less than n; otherwise the calling function should
      split the array m into suitable chunks. */
{
   mpz_t **tree;
   int *width;
   int levels;
   int i, j;

   /* Compute bottom-up a subproduct tree with m on the leaves. */
   for (i = no_m, levels = 1; i > 1; i = (i+1) / 2, levels++);
   tree = (mpz_t **) malloc (levels * sizeof (mpz_t *));
   width = (int *) malloc (levels * sizeof (int));
   width [0] = no_m;
   tree [0] = (mpz_t *) malloc (no_m * sizeof (mpz_t));
   for (j = 0; j < no_m; j++)
      mpz_init_set (tree [0][j], m [j]);
   for (i = 1; i < levels; i++) {
      width [i] = (width [i-1] + 1) / 2;
      tree [i] = (mpz_t *) malloc (width [i] * sizeof (mpz_t));
      for (j = 0; j < width [i-1] / 2; j++) {
         mpz_init (tree [i][j]);
         mpz_mul (tree [i][j], tree [i-1][2*j], tree [i-1][2*j+1]);
      }
      if (width [i-1] % 2 != 0) {
         mpz_init (tree [i][j]);
         mpz_set (tree [i][j], tree [i-1][2*j]);
      }
   }

   /* Replace top-down the tree entries by n modulo the entry. */
   mpz_mod (tree [levels-1][0], n, tree [levels-1][0]);
   for (i = levels - 2; i >= 0; i--) {
      for (j = 0; j < (width [i] / 2) * 2; j++)
         mpz_mod (tree [i][j], tree [i+1][j/2], tree [i][j]);
      if (width [i] % 2 != 0)
         mpz_set (tree [i][j], tree [i+1][j/2]);
   }

   /* Copy the leaves into the result. */
   for (j = 0; j < no_m; j++)
      mpz_set (mod [j], tree [0][j]);

   /* Clear the tree. */
   for (i = 0; i < levels; i++) {
      for (j = 0; j < width [i]; j++)
         mpz_clear (tree [i][j]);
      free (tree [i]);
   }
   free (tree);
   free (width);
}

/*****************************************************************************/

static void trial_div_batch (mpz_t *l, mpz_t *n, int no_n,
   mpz_srcptr primorialB)
   /* primorialB is supposed to be the product of the primes up to the
      smoothness bound B.
      The function removes all occurrences of the primes up to B from the
      no_n entries in n and stores the results in l, which needs to have
      the correct size and all entries of which need to be initialised. */
{
   mpz_t *mod;
   int i;

   mod = (mpz_t *) malloc (no_n * sizeof (mpz_t));
   for (i = 0; i < no_n; i++)
      mpz_init (mod [i]);

   /* Compute in mod [i] the gcd of n [i] and primorialB. */
   mpz_tree_mod (mod, primorialB, n, no_n);
   for (i = 0; i < no_n; i++)
      mpz_gcd (mod [i], n [i], mod [i]);
   /* Remove the gcd from n [i] and recompute it until all primes
      are removed. */
   for (i = 0; i < no_n; i++) {
      mpz_set (l [i], n [i]);
      while (mpz_cmp_ui (mod [i], 1ul) != 0) {
         mpz_divexact (l [i], l [i], mod [i]);
         mpz_gcd (mod [i], l [i], mod [i]);
      }
   }

   for (i = 0; i < no_n; i++)
      mpz_clear (mod [i]);
   free (mod);
}

/*****************************************************************************/

static int curve_cardinalities (mpz_t *n, mpz_srcptr N, mpz_srcptr root,
   int_cl_t d)
   /* Given N prime, a discriminant d composed of split primes modulo N,
      and a square root of d modulo N in root, the function computes the
      array of 0 (in the case that N is not a norm in Q(\sqrt d), which
      happens with probability 1 - 1 / (h/g)), 2, 4 or 6 (depending on
      the number of twists) possible cardinalities of elliptic curves
      modulo N with CM by d. The cardinalities are stored in n, the entries
      of which need to be initialised, and their number is returned. */
{
   int res;
   mpz_t t, v;
   mpz_ptr V;
   int twists;
   bool cornacchia;

   mpz_init (t);
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

   cm_timer_continue (cm_timer4);
   cm_counter4++;
   cornacchia = cm_nt_mpz_cornacchia (t, V, N, root, d);
   cm_timer_stop (cm_timer4);
   if (cornacchia) {
      res = twists;
      /* Compute the cardinalities of all the twists. */
      mpz_add_ui (n [0], N, 1);
      mpz_add (n [1], n [0], t);
      if (d == -3) {
         /* The sextic twists have trace (\pm t \pm 3*v)/2. */
         mpz_mul_ui (v, v, 3);
         mpz_sub (v, v, t);
         mpz_divexact_ui (v, v, 2);
         mpz_add (n [2], n [0], v);
         mpz_sub (n [3], n [0], v);
         mpz_add (v, v, t);
         mpz_add (n [4], n [0], v);
         mpz_sub (n [5], n [0], v);
      }
      else if (d == -4) {
         /* The quartic twists have trace \pm 2*v. */
         mpz_mul_2exp (v, v, 1);
         mpz_sub (n [2], n [0], v);
         mpz_add (n [3], n [0], v);
      }
      mpz_sub (n [0], n [0], t);
   }
   else
      res = 0;

   if (d == -3 || d == -4)
      mpz_clear (v);
   mpz_clear (t);

   return res;
}

/*****************************************************************************/

static mpz_t* compute_cardinalities (int *no_card, int_cl_t **card_d,
   mpz_srcptr N, int no_d, mpz_t *root, int_cl_t *d)
   /* Given a prime N, a list of no_d fastECPP discriminants in d and an
      array of their square roots modulo N in root, compute and return the
      array of possible CM cardinalities. The number of cardinalities is
      returned via no_card. The newly allocated array card_d contains for
      each cardinality the associated discriminant. */
{
   mpz_t *res;
   mpz_t **card;
   int *twists;
   int i, j, k;

   /* For each discriminant, compute the potential cardinalities in
      separate memory locations. */
   twists = (int *) malloc (no_d * sizeof (int));
   card = (mpz_t **) malloc (no_d * sizeof (mpz_t *));
   for (i = 0; i < no_d; i++) {
      card [i] = (mpz_t *) malloc (6 * sizeof (mpz_t));
      for (j = 0; j < 6; j++)
         mpz_init (card [i][j]);
   }

   for (i = 0; i < no_d; i++)
      twists [i] = curve_cardinalities (card [i], N, root [i], d [i]);

   /* Count the number of obtained cardinalities. */
   *no_card = 0;
   for (i = 0; i < no_d; i++)
      *no_card += twists [i];

   /* Copy the results. */
   res = (mpz_t *) malloc (*no_card * sizeof (mpz_t));
   *card_d = (int_cl_t *) malloc (*no_card * sizeof (int_cl_t));
   k = 0;
   for (i = 0; i < no_d; i++)
      for (j = 0; j < twists [i]; j++) {
         (*card_d) [k] = d [i];
         mpz_init_set (res [k], card [i][j]);
         k++;
      }

   for (i = 0; i < no_d; i++) {
      for (j = 0; j < 6; j++)
         mpz_clear (card [i][j]);
      free (card [i]);
   }
   free (card);
   free (twists);

   return res;
}

/*****************************************************************************/

static int_cl_t contains_ecpp_discriminant (mpz_ptr n, mpz_ptr l,
   mpz_srcptr N, const unsigned int delta, mpz_srcptr primorialB,
   int_cl_t *d, int no_d, mpz_t *root)
   /* The function goes through the no_d discriminants in d and tests
      whether one of them is a suitable discriminant to perform one step in
      the ECPP downrun from the (probable) prime N>=787. If one is found,
      the discriminant with the best gain is returned, and n becomes the
      cardinality of the elliptic curve and l its largest prime factor;
      otherwise 0 is returned and n and l are unchanged.
      delta >= 1 is the minimum number of bits to be gained in this step.
      primorialB is related to the trial division bound and is passed on
      to trial_div.
      root is an array of square roots of the d modulo N. */

{
   int_cl_t res;
   int no_card;
   mpz_t *card;
   int_cl_t *card_d;
   mpz_t *co;
   size_t size_co, size_N, size_opt;
   int i;

   res = 0;
   card = compute_cardinalities (&no_card, &card_d, N, no_d, root, d);

   if (no_card > 0) {
      co = (mpz_t *) malloc (no_card * sizeof (mpz_t));
      for (i = 0; i < no_card; i++)
         mpz_init (co [i]);
      cm_timer_continue (cm_timer2);
      cm_counter2++;
      trial_div_batch (co, card, no_card, primorialB);
      cm_timer_stop (cm_timer2);
      /* Look for a suitable cardinality with a point of smallest
         prime order. */
      size_N = mpz_sizeinbase (N, 2);
      size_opt = size_N;
      for (i = 0; size_opt == size_N && i < no_card; i++) {
         /* We need to check whether co > (N^1/4 + 1)^2.
            Let N have e bits, that is, 2^(e-1) <= N < 2^e,
            and co have f bits, that is, 2^(f-1) <= co < 2^f. Then
            (N^(1/4) + 1)^2 = N^(1/2) * (1+1/N^(1/4))^2
                            < 2^(e/2) * sqrt (2) for N >= 781
            So it is sufficient to check that
            f-1 >= (e+1)/2, or equivalently f >= floor (e/2) + 2.
            We also want to gain at least one bit; otherwise we might
            even increase the prime a little bit. */
         size_co = mpz_sizeinbase (co [i], 2);
         if (size_co < size_opt
             && size_co <= size_N - delta
             && size_co >= size_N / 2 + 2) {
            cm_counter3++;
            cm_timer_continue (cm_timer3);
            if (cm_nt_is_prime (co [i])) {
               res = card_d [i];
               size_opt = size_co;
               mpz_set (n, card [i]);
               mpz_set (l, co [i]);
            }
            cm_timer_stop (cm_timer3);
         }
      }
      for (i = 0; i < no_card; i++)
         mpz_clear (co [i]);
      free (co);
   }

   for (i = 0; i < no_card; i++)
      mpz_clear (card [i]);
   free (card);
   free (card_d);

   return res;
}

/*****************************************************************************/

static void root_of_d (mpz_t *Droot, int_cl_t *d, int no_d, mpz_srcptr N,
   long int *qstar, int no_qstar, mpz_t *root)
   /* Given an array of no_d discriminants d which factor over the array
      of no_qstar "signed primes" qstar, and given the roots of qstar
      modulo N in root, compute and return in Droot the roots of the d.
      Droots needs to contain enough space with all entries initialised.
      By trial dividing d over qstar, it is enough to multiply the
      corresponding roots together; however, the "even primes" need special
      care, in particular for the sign of 8. */
{
   int i, j;
   int_cl_t disc;

   for (i = 0; i < no_d; i++) {
      disc = d [i];
      mpz_set_ui (Droot [i], 1);
      for (j = 0; j < no_qstar; j++)
         if (qstar [j] % 2 != 0 && disc % qstar [j] == 0) {
            mpz_mul (Droot [i], Droot [i], root [j]);
            mpz_mod (Droot [i], Droot [i], N);
            disc /= qstar [j];
         }
      if (disc != 1) {
         for (j = 0; disc != qstar [j]; j++);
         mpz_mul (Droot [i], Droot [i], root [j]);
         mpz_mod (Droot [i], Droot [i], N);
      }
   }
}

/*****************************************************************************/

static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const unsigned int delta, mpz_srcptr primorialB,
   bool verbose, bool debug)
   /* Given a (probable) prime N>=787, return a suitable CM discriminant
      and return the cardinality of an associated elliptic curve in n and
      its largest prime factor in l.
      Dmax, hmaxprime and h are passed through to compute_discriminants.
      delta >= 1 is passed through as the minimum number of bits to be
      gained in this step.
      primorialB is passed through to trial division. */
{
   const int max_factors = 4;
   const int batch = 100;
   int no_qstar_old, no_qstar_new, no_qstar;
   long int *qstar;
   long int q;
   mpz_t *root;
   int i, j;
   int_cl_t d;
   int_cl_t *dlist;
   int no_d;
   mpz_t *Droot;

   d = 0;
   no_qstar = 0;
   q = 0;
   qstar = (long int *) malloc (0);
   root = (mpz_t *) malloc (0);
   Droot = (mpz_t *) malloc (batch * sizeof (mpz_t));
   for (i = 0; i < batch; i++)
      mpz_init (Droot [i]);

   while (d == 0) {
      /* Extend the prime and square root list. */
      cm_timer_continue (cm_timer1);
      no_qstar_old = no_qstar;
      if (no_qstar_old <= 10)
         no_qstar_new = 10;
      else
         no_qstar_new = 20;
      no_qstar = no_qstar_old + no_qstar_new;
      qstar = (long int *) realloc (qstar, no_qstar * sizeof (long int));
      root = (mpz_t *) realloc (root, no_qstar * sizeof (mpz_t));
      for (i = no_qstar_old; i < no_qstar; i++)
         mpz_init (root [i]);
      compute_qstar (qstar + no_qstar_old, root + no_qstar_old, N, &q,
         no_qstar_new);
      cm_timer_stop (cm_timer1);

      /* Precompute a list of potential discriminants for fastECPP. */
      dlist = compute_sorted_discriminants (&no_d, qstar, no_qstar_old,
         no_qstar_new, max_factors, Dmax, hmaxprime, h);
      if (verbose && debug)
         printf ("%6i discriminants for %3i primes\n", no_d, no_qstar);

      /* Go through the list, treating batch discriminants at a time. */
      for (i = 0; d == 0 && i < (no_d + batch - 1) / batch; i++) {
         j = no_d - i * batch;
         if (j > batch)
            j = batch;
         root_of_d (Droot, dlist + i * batch, j, N, qstar, no_qstar, root);
         d = contains_ecpp_discriminant (n, l, N, delta, primorialB,
            dlist + i * batch, j, Droot);
      }

      /* Free the discriminant list. */
      free (dlist);
   }

   for (i = 0; i < batch; i++)
      mpz_clear (Droot [i]);
   free (Droot);
   for (i = 0; i < no_qstar; i++)
      mpz_clear (root [i]);
   free (root);
   free (qstar);

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
   const size_t L = mpz_sizeinbase (p, 2);
   const unsigned long int B = (L >> 5) * (L >> 5) * (L >> 5);
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
   cm_timer_t clock;
   int i;

   /* Precompute class numbers. */
   cm_timer_start (clock);
   Dmax = ((L * L) >> 4) << 2;
   h = (uint_cl_t *) malloc ((Dmax / 2) * sizeof (uint_cl_t));
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
      printf ("-- Time for primorial of B=%lu: %5.1f\n", B,
         cm_timer_get (clock));

   mpz_init_set (N, p);
   *depth = 0;
   c = (mpz_t**) malloc (*depth);
   cm_timer_reset (cm_timer1);
   cm_timer_reset (cm_timer2);
   cm_timer_reset (cm_timer3);
   cm_timer_reset (cm_timer4);
   cm_counter1 = 0;
   cm_counter2 = 0;
   cm_counter3 = 0;
   cm_counter4 = 0;
   while (mpz_sizeinbase (N, 2) > 64) {
      c = (mpz_t**) realloc (c, (*depth + 1) * sizeof (mpz_t *));
      c [*depth] = (mpz_t *) malloc (4 * sizeof (mpz_t));
      for (i = 0; i < 4; i++)
         mpz_init (c [*depth][i]);
      mpz_set (c [*depth][0], N);
      cm_timer_start (clock);
      if (verbose)
         printf ("-- Size %4li bits\n", mpz_sizeinbase (N, 2));
      d = find_ecpp_discriminant (c [*depth][2], c [*depth][3], N, Dmax,
         hmaxprime, h, delta, primorialB, verbose, debug);
      cm_timer_stop (clock);
      if (verbose) {
         printf ("   Time for discriminant %8"PRIicl": %5.1f\n",
            d, cm_timer_get (clock));
         if (debug) {
            printf ("   largest prime of d: %"PRIucl"\n",
                    cm_nt_largest_factor (-d));
            printf ("   largest prime of h: %"PRIucl"\n",
                    cm_nt_largest_factor (h [(-d) / 2 - 1]));
            printf ("%6i qstar:      %.1f\n",
                  cm_counter1, cm_timer_get (cm_timer1));
            printf ("%6i Trial div:  %.1f\n", cm_counter2,
                  cm_timer_get (cm_timer2));
            printf ("%6i is_prime:   %.1f\n", cm_counter3,
                  cm_timer_get (cm_timer3));
            printf ("%6i Cornacchia: %.1f\n", cm_counter4,
                  cm_timer_get (cm_timer4));
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

static void ecpp_param_init (cm_param_ptr param, uint_cl_t d)
   /* Given the discriminant d, find the best parameter combination in the
      ECPP setting. To minimise the number of polynomial factorisations, we
      only use invariants that do not require modular polynomials, and
      additionally a small number of invariants with low-degree modular
      polynomials. Among these, we use the one with the best height
      factor. */
{
   cm_param_t new_param;
   double hf, new_hf;

   /* First test the non-parametric invariants in the good order. */
   if (   !cm_param_init (param, d, CM_INVARIANT_WEBER,
              CM_SUBFIELD_PREFERRED, 0, false)
       && !(((d - 1) % 8 == 0
            && cm_param_init (param, 4*d, CM_INVARIANT_WEBER,
                  0, CM_SUBFIELD_PREFERRED, false)))
       && !cm_param_init (param, d, CM_INVARIANT_GAMMA2,
              0, CM_SUBFIELD_PREFERRED, false)
       && !cm_param_init (param, d, CM_INVARIANT_GAMMA3,
              0, CM_SUBFIELD_PREFERRED, false))
       cm_param_init (param, d, CM_INVARIANT_J,
              0, CM_SUBFIELD_PREFERRED, false);
   hf = cm_class_height_factor (param);

   /* Atkin invariants have excellent factors between 24 and 36, but
      the Hecke operators are slow to compute. So do not use them. */

   /* Simple eta uses the best of the w_n^e with n one of
      3, 5, 7, 13, 4, 9 or 25, the values for which the modular
      polynomial has genus 0. */
   if (cm_param_init (new_param, d, CM_INVARIANT_SIMPLEETA,
          0, CM_SUBFIELD_PREFERRED, false)) {
      new_hf = cm_class_height_factor (new_param);
      if (new_hf > hf) {
         param [0] = new_param [0];
         hf = new_hf;
      }
   }

   /* For double eta quotients, we limit the search to a degree
      of 2 in j. */
   if (cm_param_init (new_param, d, CM_INVARIANT_DOUBLEETA,
          -1, CM_SUBFIELD_PREFERRED, false)) {
      new_hf = cm_class_height_factor (new_param);
      if (new_hf > hf) {
         param [0] = new_param [0];
         hf = new_hf;
      }
   }

   /* Multiple eta quotients slow the code down; this is probably due
      to the need for factoring the modular polynomials of degree 4. */
}

/*****************************************************************************/

static void ecpp2_one_step (mpz_t *cert2, mpz_t *cert1,
   const char* modpoldir, bool tower, bool verbose, bool debug)
   /* The function takes the same parameters as cm_ecpp2, except that cert1
      contains only one entry of the first ECPP step, and that only one
      step of the certificate is output in cert2, which needs be to a
      pre-allocated and initialised array with 6 entries. This function
      has been separated to help with parallelisation. */
{
   mpz_srcptr p, n, l;
   int_cl_t d;
   mpz_t t, co, a, b, x, y;
   cm_param_t param;
   cm_class_t c;
   cm_timer_t clock;

   mpz_init (t);
   mpz_init (co);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);

   cm_timer_start (clock);
   p = cert1 [0];
   d = mpz_get_si (cert1 [1]);
   n = cert1 [2];
   l = cert1 [3];

   mpz_add_ui (t, p, 1);
   mpz_sub (t, t, n);
   mpz_divexact (co, n, l);

   cm_timer_continue (cm_timer3);
   ecpp_param_init (param, d);

   /* Compute one of the class field tower or the class polynomial. */
   cm_class_init (c, param, false);
   cm_class_compute (c, param, !tower, tower, false);
   cm_timer_stop (cm_timer3);

   cm_curve_and_point (a, b, x, y, param, c, p, l, co,
      modpoldir, false, false);
   cm_class_clear (c);
   cm_timer_stop (clock);

   if (verbose) {
      printf ("-- Time for discriminant %8"PRIicl" with invariant %c "
         "and parameters %s for %4li bits: %5.1f\n",
         d, param->invariant, param->str, mpz_sizeinbase (p, 2),
         cm_timer_get (clock));
      if (debug) {
         printf ("   CM:    %5.1f\n", cm_timer_get (cm_timer3));
         printf ("   roots: %5.1f\n", cm_timer_get (cm_timer1));
         printf ("   point: %5.1f\n", cm_timer_get (cm_timer2));
      }
   }

   mpz_set (cert2 [0], p);
   mpz_set (cert2 [1], t);
   mpz_set (cert2 [2], co);
   mpz_set (cert2 [3], a);
   mpz_set (cert2 [4], x);
   mpz_set (cert2 [5], y);

   mpz_clear (t);
   mpz_clear (co);
   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (x);
   mpz_clear (y);
}

/*****************************************************************************/

static void cm_ecpp2 (mpz_t **cert2, mpz_t **cert1, int depth,
   const char* modpoldir, bool tower, bool verbose, bool debug)
   /* Given the result of the ECPP down-run in cert1, an array of length
      depth as computed by cm_ecpp1 or cm_pari_ecpp1, execute the second
      step of the ECPP algorithm and compute the certificate proper in
      cert2, which needs to be pre-allocated as an array of length depth
      with each entry an array of length 6 to contain p, t, co, a, x and y:
      the prime p to be certified, the trace t of the elliptic curve, the
      co-factor such that p+1-t = co*l with l the next prime, the curve
      parameter a, and the coordinates (x,y) of a point of prime order l
      on the curve; the curve parameter b is implicit.
      modpoldir gives the directory where modular polynomials are stored;
      it is passed through to the function computing a curve from a root
      of the class polynomial.
      tower indicates whether a class field tower decomposition is used
      instead of only the class polynomial.
      verbose indicates whether intermediate computations output progress
      information.
      debug indicates whether additional developer information (mainly
      timings and counters for tuning) is output; this is done only in
      the case that verbose is set as well. */

{
   int i;

   cm_timer_reset (cm_timer1);
   cm_timer_reset (cm_timer2);
   cm_timer_reset (cm_timer3);

   for (i = 0; i < depth; i++)
      ecpp2_one_step (cert2 [i], cert1 [i], modpoldir, tower,
         verbose, debug);
}

/*****************************************************************************/

bool cm_ecpp (mpz_srcptr N, const char* modpoldir, bool pari, bool tower,
   bool print, bool check, bool verbose, bool debug)
   /* Assuming that N is a (probable) prime, compute an ECPP certificate.
      modpoldir gives the directory where modular polynomials are stored;
      it is passed through to the function computing a curve from a root
      of the class polynomial.
      pari indicates whether the first, downrun step of PARI is used
      instead of the built-in function.
      tower indicates whether a class field tower decomposition is used
      instead of only the class polynomial.
      print indicates whether the result is printed.
      check indicates whether the certificate should be checked.
      If yes, the return value of the function is the result of the check;
      otherwise the return value is true.
      verbose indicates whether intermediate computations output progress
      information.
      debug indicates whether additional developer information (mainly
      timings and counters for tuning) is output; this is done only in
      the case that verbose is set as well. */

{
   bool res = true;
   int depth;
   mpz_t **cert1, **cert2;
   int i, j;
   cm_timer_t clock, clock2;

   cm_timer_start (clock);
   cm_timer_start (clock2);
   if (pari)
      cert1 = cm_pari_ecpp1 (&depth, N);
   else
      cert1 = cm_ecpp1 (&depth, N, verbose, debug);
   cm_timer_stop (clock2);
   if (verbose)
      printf ("--- Time for first ECPP step, depth %i:  %.1f\n", depth,
         cm_timer_get (clock2));

   cm_timer_start (clock2);
   cert2 = (mpz_t **) malloc (depth * sizeof (mpz_t *));
   for (i = 0; i < depth; i++) {
      cert2 [i] = (mpz_t *) malloc (6 * sizeof (mpz_t));
      for (j = 0; j < 6; j++)
         mpz_init (cert2 [i][j]);
   }
   cm_ecpp2 (cert2, cert1, depth, modpoldir, tower, verbose, debug);
   if (print) {
      printf ("c = [");
      for (i = 0; i < depth; i++) {
         printf ("[");
         for (j = 0; j < 4; j++) {
            mpz_out_str (stdout, 10, cert2 [i][j]);
            printf (", ");
         }
         printf ("[");
         mpz_out_str (stdout, 10, cert2 [i][4]);
         printf (", ");
         mpz_out_str (stdout, 10, cert2 [i][5]);
         printf ("]]");
         if (i != depth - 1)
            printf (", ");
      }
      printf ("];\n");
   }

   cm_timer_stop (clock2);
   cm_timer_stop (clock);
   if (verbose) {
      printf ("--- Time for second ECPP step: %.1f\n",
         cm_timer_get (clock2));
      printf ("--- Total time for ECPP:       %.1f\n",
         cm_timer_get (clock));
   }

   if (check) {
      cm_timer_start (clock);
      res = cm_pari_ecpp_check (cert2, depth);
      cm_timer_stop (clock);
      if (verbose)
         printf ("--- Time for ECPP check (%s): %.1f\n",
            (res ? "true" : "false"), cm_timer_get (clock));
   }

   for (i = 0; i < depth; i++) {
      for (j = 0; j < 4; j++)
         mpz_clear (cert1 [i][j]);
      for (j = 0; j < 6; j++)
         mpz_clear (cert2 [i][j]);
      free (cert1 [i]);
      free (cert2 [i]);
   }
   free (cert1);
   free (cert2);

   return res;
}

/*****************************************************************************/
/*****************************************************************************/
