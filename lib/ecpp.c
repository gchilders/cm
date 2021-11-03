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

static void compute_h (uint_cl_t *h, uint_cl_t Dmax, cm_stat_t stat);
static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   long int *q, int no, cm_stat_t stat);
static int_cl_t** compute_signed_discriminants (int *no_d, long int *qstar,
   int no_qstar, uint_cl_t Dmax, int sign);
static int_cl_t** compute_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h);
static int_cl_t* compute_sorted_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const double prob, bool debug);
static int disc_cmp (const void* d1, const void* d2);
static void mpz_tree_mod (mpz_t *mod, mpz_srcptr n, mpz_t *m, int no_m);
static void trial_div_batch (mpz_t *l, mpz_t *n, int no_n,
   mpz_srcptr primorialB);
static mpz_t* compute_cardinalities (int *no_card, int_cl_t **card_d,
   mpz_srcptr N, int no_d, mpz_t *root, int_cl_t *d, cm_stat_ptr stat);
static int card_cmp (const void* c1, const void* c2);
static int_cl_t contains_ecpp_discriminant (mpz_ptr n, mpz_ptr l,
   mpz_srcptr N, mpz_t *card, mpz_t *l_list, int_cl_t *d, int no_card,
   const unsigned int delta, cm_stat_t stat);
static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const unsigned int delta, mpz_srcptr primorialB,
   bool debug, cm_stat_t stat);
static void ecpp_param_init (cm_param_ptr param, uint_cl_t d);
static mpz_t** cm_ecpp1 (int *depth, mpz_srcptr p, bool verbose,
   bool debug, cm_stat_ptr stat);
static void cm_ecpp2 (mpz_t **cert2, mpz_t **cert1, int depth,
   const char* modpoldir, bool tower, bool verbose, bool debug,
   cm_stat_ptr stat);

/*****************************************************************************/

void cm_ecpp_compute_h_chunk (uint_cl_t *h, uint_cl_t Dmin, uint_cl_t Dmax)
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
      printf ("***** Error: cm_ecpp_compute_h_chunk called with "
         "parameters not divisible by 4.\n");
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

static void compute_h (uint_cl_t *h, uint_cl_t Dmax, cm_stat_t stat)
   /* The function behaves as compute_h_chunk (h, 0, Dmax), but computing
      the class numbers in ranges of 100000 (that is, 50000 discriminants
      at a time). Experimentally, this optimises the running time on my
      laptop for Dmax=10^7. The optimal value probably depends on the cache
      size. */
{
   const uint_cl_t size = 100000;
   uint_cl_t last;
   int chunks;
#ifdef WITH_MPI
   MPI_Status status;
   int sent, received, rank, job;
   double t;
   cm_stat_t stat_worker;
#else
   int i;
#endif

   if (Dmax % 4 != 0) {
      printf ("***** Error: compute_h called with parameter not "
              "divisible by 4.\n");
      exit (1);
   }

   chunks = (Dmax + size - 1) / size;
#ifndef WITH_MPI
   cm_timer_start (stat->timer [5]);
   for (i = 0; i < chunks; i++) {
      last = (i + 1) * size;
      if (last > Dmax)
         last = Dmax;
      cm_ecpp_compute_h_chunk (h + i * size / 2, i * size, last);
   }
   cm_timer_stop (stat->timer [5]);
#else
   sent = 0;
   received = 0;
   t = cm_timer_get (stat->timer [5]);
   cm_timer_continue (stat->timer [5]);
   while (received < chunks) {
      if (sent < chunks && (rank = cm_mpi_queue_pop ()) != -1) {
         last = (sent + 1) * size;
         if (last > Dmax)
            last = Dmax;
         cm_mpi_submit_h_chunk (rank, sent, sent * size, last);
         sent++;
      }
      else {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         cm_mpi_get_h_chunk (h + job * size / 2, rank, stat_worker);
         t += cm_timer_get (stat_worker->timer [5]);
         cm_mpi_queue_push (rank);
         received++;
      }
   }
   cm_timer_stop (stat->timer [5]);
   stat->timer [5]->elapsed = t;
#endif
}

/*****************************************************************************/

static void compute_qstar (long int *qstar, mpz_t *root, mpz_srcptr p,
   long int *q, int no, cm_stat_t stat)
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
#ifdef WITH_MPI
   MPI_Status status;
   int sent, received, rank, job;
   double t;
   cm_stat_t stat_worker;
#endif

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
#ifndef WITH_MPI
   cm_timer_continue (stat->timer [0]);
   for (i = 0; i < no; i++)
      cm_nt_mpz_tonelli_si_with_generator (root [i], qstar [i], p, e, r, z);
   cm_timer_stop (stat->timer [0]);
#else
   sent = 0;
   received = 0;
   t = cm_timer_get (stat->timer [0]);
      /* Memorise CPU time to avoid confusion with server CPU time. */
   cm_timer_continue (stat->timer [0]);
   while (received < no) {
      if (sent < no && (rank = cm_mpi_queue_pop ()) != -1) {
         cm_mpi_submit_tonelli (rank, sent, qstar [sent], p, e, r, z);
         sent++;
      }
      else {
         /* There is either nothing to send, or all workers are busy.
            Wait for a result. */
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         cm_mpi_get_tonelli (root [job], rank, stat_worker);
         t += cm_timer_get (stat_worker->timer [0]);
         cm_mpi_queue_push (rank);
         received++;
      }
   }
   cm_timer_stop (stat->timer [0]);
   stat->timer [0]->elapsed = t;
      /* Restore CPU time containing all the workers, while keeping
         the wallclock time as measured by the server. */
#endif

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

static int_cl_t* compute_sorted_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const double prob, bool debug)
   /* The function takes the same parameters as compute_discriminants (and
      most of them are just passed through). But instead of an unsorted
      double array, it returns a simple array of discriminants sorted
      according to the usefulness in ECPP.
      The argument prob serves for debugging purposes; it gives the
      probability, computed from the size of the number to be proved prime
      and the trial division bound, that a curve cardinality will lead to
      success in this step. */
{
   int_cl_t **dlist;
   int_cl_t *d;
   double exp_card;
      /* the expected number of curve cardinalities; the sum over
         #twists * (g/h) */
   int i;

   dlist = compute_discriminants (no_d, qstar, no_qstar_old, no_qstar_new,
      max_factors, Dmax, hmaxprime, h);
   qsort (dlist, *no_d, sizeof (int_cl_t *), disc_cmp);

   exp_card = 0;
   d = (int_cl_t *) malloc (*no_d * sizeof (int_cl_t));
   for (i = 0; i < *no_d; i++) {
      d [i] = dlist [i][0];
      if (d [i] == -3)
         exp_card += 6.0;
      else if (d [i] == -4)
         exp_card += 4.0;
      else
         exp_card += 2.0 / dlist [i][2];
      free (dlist [i]);
   }
   free (dlist);

   if (debug)
      printf ("   expected number of primes: %.1f\n", exp_card * prob);

   return d;
}

/*****************************************************************************/

void cm_ecpp_sqrt_d (mpz_t *Droot, int_cl_t *d, int no_d, mpz_srcptr N,
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

int cm_ecpp_curve_cardinalities (mpz_t *n, mpz_srcptr N, mpz_srcptr root,
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

   cornacchia = cm_pari_cornacchia (t, V, N, root, d);
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
   mpz_srcptr N, int no_d, mpz_t *root, int_cl_t *d, cm_stat_ptr stat)
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
#ifdef WITH_MPI
   MPI_Status status;
   int sent, received, rank, job;
   double t;
   cm_stat_t stat_worker;
#endif


   /* For each discriminant, compute the potential cardinalities in
      separate memory locations. */
   twists = (int *) malloc (no_d * sizeof (int));
   card = (mpz_t **) malloc (no_d * sizeof (mpz_t *));
   for (i = 0; i < no_d; i++) {
      card [i] = (mpz_t *) malloc (6 * sizeof (mpz_t));
      for (j = 0; j < 6; j++)
         mpz_init (card [i][j]);
   }

#ifndef WITH_MPI
   cm_timer_continue (stat->timer [3]);
   for (i = 0; i < no_d; i++)
      twists [i] = cm_ecpp_curve_cardinalities (card [i], N, root [i],
         d [i]);
   cm_timer_stop (stat->timer [3]);
#else
   sent = 0;
   received = 0;
   t = cm_timer_get (stat->timer [3]);
   cm_timer_continue (stat->timer [3]);
   while (received < no_d) {
      if (sent < no_d && (rank = cm_mpi_queue_pop ()) != -1) {
         cm_mpi_submit_curve_cardinalities (rank, sent, N, root [sent],
            d [sent]);
         sent++;
      }
      else {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         twists [job] = cm_mpi_get_curve_cardinalities (card [job], rank,
            stat_worker);
         t += cm_timer_get (stat_worker->timer [3]);
         cm_mpi_queue_push (rank);
         received++;
      }
   }
   cm_timer_stop (stat->timer [3]);
   stat->timer [3]->elapsed = t;
#endif
   stat->counter [3] += no_d;

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

static int card_cmp (const void* c1, const void* c2)
{
   mpz_t *C1, *C2;

   C1 = (*((mpz_t **) c1));
   C2 = (*((mpz_t **) c2));

   /* Sort by the first component. */
   return mpz_cmp (C1 [0], C2 [0]);
}

/*****************************************************************************/

static int_cl_t contains_ecpp_discriminant (mpz_ptr n, mpz_ptr l,
   mpz_srcptr N, mpz_t *card, mpz_t *l_list, int_cl_t *d, int no_card,
   const unsigned int delta, cm_stat_t stat)
   /* For the no_card discriminants in d, card is supposed to contain
      corresponding curve cardinalities and l_list their non-smooth parts.
      The function tests whether one of them is suitable to perform one
      step in the ECPP downrun from the (probable) prime N>=787. If one is
      found, the corresponding discriminant from d is returned, and n
      becomes the cardinality of the elliptic curve and l its largest prime
      factor; otherwise 0 is returned and n and l are unchanged.
      delta >= 1 is the minimum number of bits to be gained in this
      step. */

{
   int_cl_t res;
   size_t size_l, size_N;
   int no, batch;
   int *index;
   mpz_t **c;
   int i, j, max_i, max_j;
#ifdef WITH_MPI
   MPI_Status status;
   int sent, received, size, rank, job;
   double t;
   cm_stat_t stat_worker;
#endif

   res = 0;
   size_N = mpz_sizeinbase (N, 2);

   /* Filter out suitable cardinalities and copy them to a new array;
      each entry is a 2-dimensional array containing the potential prime
      factor and its original index, so that the cardinality and the
      discriminant can be found again after sorting. */
   c = (mpz_t **) malloc (no_card * sizeof (mpz_t *));
   no = 0;
   for (i = 0; i < no_card; i++) {
      /* We need to check whether l > (N^1/4 + 1)^2.
         Let N have e bits, that is, 2^(e-1) <= N < 2^e,
         and l have f bits, that is, 2^(f-1) <= l < 2^f. Then
         (N^(1/4) + 1)^2 = N^(1/2) * (1+1/N^(1/4))^2
         < 2^(e/2) * sqrt (2) for N >= 781.
         So it is sufficient to check that
         f-1 >= (e+1)/2, or equivalently f >= floor (e/2) + 2. */
      size_l = mpz_sizeinbase (l_list [i], 2);
      if (   size_l <= size_N - delta
          && size_l >= size_N / 2 + 2) {
         c [no] = (mpz_t *) malloc (2 * sizeof (mpz_t));
         mpz_init_set (c [no][0], l_list [i]);
         mpz_init_set_ui (c [no][1], i);
         no++;
      }
   }
   c = (mpz_t **) realloc (c, no * sizeof (mpz_t *));

   /* Sort by increasing non-smooth part. */
   qsort (c, no, sizeof (mpz_t *), card_cmp);

   /* Go through the array in batches of size batch and look for a prime
      non-smooth part. Stop and remember the first occurrence when it is
      found; this is the smallest possible prime in the list.
      Using batches is only meaningful in the case of parallelisation. */
#ifndef WITH_MPI
   batch = 1;
#else
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   batch = size - 1;
#endif
   index = (int *) malloc (batch * sizeof (int));
   max_i = (no + batch - 1) / batch;
   for (i = 0; res == 0 && i < max_i; i++) {
      max_j = (i + 1) * batch;
      if (max_j > no)
         max_j = no;
      for (j = 0; j < batch; j++)
         index [j] = -1;
#ifndef WITH_MPI
      cm_timer_continue (stat->timer [2]);
      for (j = i * batch; j < max_j; j++) {
         if (cm_nt_is_prime (c [j][0]))
            index [j - i * batch] = mpz_get_ui (c [j][1]);
      }
      cm_timer_stop (stat->timer [2]);
#else
   sent = i * batch;
   received = 0;
   t = cm_timer_get (stat->timer [2]);
   cm_timer_continue (stat->timer [2]);
   while (received < max_j - i * batch) {
      if (sent < max_j && (rank = cm_mpi_queue_pop ()) != -1) {
         cm_mpi_submit_is_prime (rank, sent, c [sent][0]);
         sent++;
      }
      else {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         if (cm_mpi_get_is_prime (rank, stat_worker))
            index [job - i * batch] = mpz_get_ui (c [job][1]);
         t += cm_timer_get (stat_worker->timer [0]);
         cm_mpi_queue_push (rank);
         received++;
      }
   }
   cm_timer_stop (stat->timer [2]);
   stat->timer [2]->elapsed = t;
#endif
      for (j = 0; res == 0 && j < batch; j++)
         if (index [j] != -1) {
            res = d [index [j]];
            mpz_set (n, card [index [j]]);
            mpz_set (l, l_list [index [j]]);
         }
      stat->counter [2] += max_j - i * batch;
   }
   free (index);

   for (i = 0; i < no; i++) {
      mpz_clear (c [i][0]);
      mpz_clear (c [i][1]);
      free (c [i]);
   }
   free (c);

   return res;
}

/*****************************************************************************/

static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, uint_cl_t *h,
   const unsigned int delta, mpz_srcptr primorialB,
   bool debug, cm_stat_t stat)
   /* Given a (probable) prime N>=787, return a suitable CM discriminant
      and return the cardinality of an associated elliptic curve in n and
      its largest prime factor in l.
      Dmax, hmaxprime and h are passed through to compute_discriminants.
      delta >= 1 is passed through as the minimum number of bits to be
      gained in this step.
      primorialB is passed through to trial division. */
{
   const int max_factors = 4;
   const int no_qstar_delta [] = { 20, 40, 50 };
      /* Number of new qstar to add, roughly optimised through removing the
         first 1000 bits of nextprime (10^(1000*i)). */
   int no_qstar_old, no_qstar_new, no_qstar;
   long int *qstar;
   long int q;
   mpz_t *root, *Droot, *card, *l_list;
   int_cl_t d;
   int_cl_t *dlist, *d_card;
   int no_d, no_card, batch, no_d_batch;
   int i, max_i;
   const double prob = 1.7811
      * log2 (mpz_sizeinbase (primorialB, 2) * M_LN2)
      / mpz_sizeinbase (N, 2);
      /* According to [FrKlMoWi04], the probability that a number N is
         a B-smooth part times a prime is exp (gamma) * log B / log N.
         Here we do not know B, but it is quite precisely the logarithm
         of its primorial. */
#ifdef WITH_MPI
   MPI_Status status;
   int sent, received, size, rank, job;
   double t;
   cm_stat_t stat_worker;
#endif

   d = 0;
   no_qstar = 0;
   q = 0;
   qstar = (long int *) malloc (0);
   root = (mpz_t *) malloc (0);

   i = mpz_sizeinbase (N, 2) / 3322;
   if (i >= (int) (sizeof (no_qstar_delta) / sizeof (no_qstar_delta [0])))
      i = sizeof (no_qstar_delta) / sizeof (no_qstar_delta [0]) - 1;
   no_qstar_new = no_qstar_delta [i];

   while (d == 0) {
      /* Extend the prime and square root list. */
      no_qstar_old = no_qstar;
      no_qstar = no_qstar_old + no_qstar_new;
      qstar = (long int *) realloc (qstar, no_qstar * sizeof (long int));
      root = (mpz_t *) realloc (root, no_qstar * sizeof (mpz_t));
      for (i = no_qstar_old; i < no_qstar; i++)
         mpz_init (root [i]);
      compute_qstar (qstar + no_qstar_old, root + no_qstar_old, N, &q,
            no_qstar_new, stat);
      stat->counter [0] += no_qstar_new;

      /* Precompute a list of potential discriminants. This takes
         virtually no time, so we do not measure it.*/
      dlist = compute_sorted_discriminants (&no_d, qstar, no_qstar_old,
            no_qstar_new, max_factors, Dmax, hmaxprime, h, prob, debug);

      /* Precompute their square roots modulo N in batches, which is useful
         for parallelisation. */
      Droot = (mpz_t *) malloc (no_d * sizeof (mpz_t));
      for (i = 0; i < no_d; i++)
         mpz_init (Droot [i]);
#ifndef WITH_MPI
      batch = no_d;
#else
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      batch = (no_d + size - 2) / (size - 1);
#endif
      max_i = (no_d + batch - 1) / batch;
#ifndef WITH_MPI
      cm_timer_continue (stat->timer [4]);
      for (i = 0; i < max_i; i++) {
         no_d_batch = no_d - i * batch;
         if (no_d_batch > batch)
            no_d_batch = batch;
         cm_ecpp_sqrt_d (Droot + i * batch, dlist + i * batch, no_d_batch,
               N, qstar, no_qstar, root);
      }
      cm_timer_stop (stat->timer [4]);
#else
      sent = 0;
      received = 0;
      t = cm_timer_get (stat->timer [4]);
      cm_timer_continue (stat->timer [4]);
      while (received < max_i) {
         if (sent < max_i && (rank = cm_mpi_queue_pop ()) != -1) {
            no_d_batch = no_d - sent * batch;
            if (no_d_batch > batch)
               no_d_batch = batch;
            cm_mpi_submit_sqrt_d (rank, sent, dlist + sent * batch,
               no_d_batch, N, qstar, no_qstar, root);
            sent++;
         }
         else {
            MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
                  MPI_COMM_WORLD, &status);
            rank = status.MPI_SOURCE;
            cm_mpi_get_sqrt_d (Droot + job * batch, rank, stat_worker);
            t += cm_timer_get (stat_worker->timer [4]);
            cm_mpi_queue_push (rank);
            received++;
         }
      }
      cm_timer_stop (stat->timer [4]);
      stat->timer [4]->elapsed = t;
#endif
      stat->counter [4] += no_d;

      /* Compute the cardinalities of the corresponding elliptic curves. */
      card = compute_cardinalities (&no_card, &d_card, N,
            no_d, Droot, dlist, stat);

      if (no_card > 0) {
         /* Remove smooth parts of cardinalities. */
         l_list = (mpz_t *) malloc (no_card * sizeof (mpz_t));
         for (i = 0; i < no_card; i++)
            mpz_init (l_list [i]);
         cm_timer_continue (stat->timer [1]);
         stat->counter [1]++;
         trial_div_batch (l_list, card, no_card, primorialB);
         cm_timer_stop (stat->timer [1]);

         d = contains_ecpp_discriminant (n, l, N, card, l_list, d_card,
               no_card, delta, stat);

         for (i = 0; i < no_card; i++)
            mpz_clear (l_list [i]);
         free (l_list);
      }

      free (dlist);
      for (i = 0; i < no_d; i++)
         mpz_clear (Droot [i]);
      free (Droot);
      for (i = 0; i < no_card; i++)
         mpz_clear (card [i]);
      free (card);
   }

   for (i = 0; i < no_qstar; i++)
      mpz_clear (root [i]);
   free (root);
   free (qstar);

   return d;
}

/*****************************************************************************/

static mpz_t** cm_ecpp1 (int *depth, mpz_srcptr p, bool verbose,
   bool debug, cm_stat_ptr stat)
   /* Compute the first step of the ECPP certificate; this is the downrun
      part with the parameters of the elliptic curves.
      The return value is a newly allocated array of depth entries, each
      of which is an array of length 4, containing in this order
      - p_i, a prime to be certified;
      - d_i, the discriminant;
      - n_i, the cardinality of the elliptic curve;
      - l_i, the prime order dividing this cardinality.
      The downrun stops as soon as the prime is less than 2^64.
      The stat parameter makes it possible to pass timing information
      to the calling function. */

{
   const size_t L = mpz_sizeinbase (p, 2);
   const unsigned long int B = (L >> 4) * (L >> 4) * (L >> 5);
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
   double t, t_old;
   int i;

   cm_stat_init (stat);
   cm_timer_start (stat->timer [7]);

   /* Precompute class numbers. */
   Dmax = ((L * L) >> 4) << 2;
   h = (uint_cl_t *) malloc ((Dmax / 2) * sizeof (uint_cl_t));
   compute_h (h, Dmax, stat);
   if (verbose)
      printf ("-- Time for class numbers up to Dmax=%"PRIucl
         ": %.1f (%.1f)\n", Dmax,
         cm_timer_get (stat->timer [5]),
         cm_timer_wc_get (stat->timer [5]));

   /* Precompute primorial for trial division. */
   cm_timer_start (stat->timer [6]);
   mpz_init (primorialB);
   mpz_primorial_ui (primorialB, B);
   cm_timer_stop (stat->timer [6]);
   if (verbose)
      printf ("-- Time for primorial of B=%lu: %5.1f\n", B,
         cm_timer_get (stat->timer [6]));

   mpz_init_set (N, p);
   *depth = 0;
   c = (mpz_t**) malloc (*depth);
   t_old = 0;
   while (mpz_sizeinbase (N, 2) > 64) {
      c = (mpz_t**) realloc (c, (*depth + 1) * sizeof (mpz_t *));
      c [*depth] = (mpz_t *) malloc (4 * sizeof (mpz_t));
      for (i = 0; i < 4; i++)
         mpz_init (c [*depth][i]);
      mpz_set (c [*depth][0], N);
      cm_timer_start (clock);
      if (verbose)
         printf ("-- Size [%i]: %4li bits\n", *depth,
            mpz_sizeinbase (N, 2));
      d = find_ecpp_discriminant (c [*depth][2], c [*depth][3], N, Dmax,
         hmaxprime, h, delta, primorialB, debug, stat);
      cm_timer_stop (clock);
      if (verbose) {
         for (i = 0, t = 0.0; i < 5; i++)
            t += cm_timer_get (stat->timer [i]);
         printf ("   Time for discriminant %8"PRIicl": %5.1f (%5.1f)\n",
            d, t - t_old, cm_timer_wc_get (clock));
         t_old = t;
         if (debug) {
            printf ("   largest prime of d: %"PRIucl"\n",
               cm_nt_largest_factor (-d));
            printf ("   largest prime of h: %"PRIucl"\n",
               cm_nt_largest_factor (h [(-d) / 2 - 1]));
            printf ("%6i qstar:      %.1f (%.1f)\n", stat->counter [0],
               cm_timer_get (stat->timer [0]),
               cm_timer_wc_get (stat->timer [0]));
            printf ("%6i sqrt:       %.1f (%.1f)\n", stat->counter [4],
               cm_timer_get (stat->timer [4]),
               cm_timer_wc_get (stat->timer [4]));
            printf ("%6i Cornacchia: %.1f (%.1f)\n", stat->counter [3],
               cm_timer_get (stat->timer [3]),
               cm_timer_wc_get (stat->timer [3]));
            printf ("%6i Trial div:  %.1f\n", stat->counter [1],
               cm_timer_get (stat->timer [1]));
            printf ("%6i is_prime:   %.1f (%.1f)\n", stat->counter [2],
                  cm_timer_get (stat->timer [2]),
                  cm_timer_wc_get (stat->timer [2]));
         }
      }
      mpz_set_si (c [*depth][1], d);
      mpz_set (N, c [*depth][3]);
      (*depth)++;
   }

   free (h);
   mpz_clear (primorialB);
   mpz_clear (N);

   cm_timer_stop (stat->timer [7]);
   for (i = 0, t = 0.0; i < 7; i++)
      t+= cm_timer_get (stat->timer [i]);
   if (verbose)
      printf ("--- Time for first ECPP step, depth %i:  %.1f (%.1f)\n",
         *depth, t, cm_timer_wc_get (stat->timer [7]));

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

void cm_ecpp_one_step2 (mpz_t *cert2, mpz_t *cert1,
   const char* modpoldir, bool tower, bool verbose, bool debug,
   cm_stat_t stat)
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

   cm_timer_continue (stat->timer [1]);
   ecpp_param_init (param, d);

   /* Compute one of the class field tower or the class polynomial. */
   cm_class_init (c, param, false);
   cm_class_compute (c, param, !tower, tower, false);
   cm_timer_stop (stat->timer [1]);

   cm_curve_and_point_stat (a, b, x, y, param, c, p, l, co,
      modpoldir, false, false, stat);
   cm_class_clear (c);
   cm_timer_stop (clock);

   if (verbose) {
      printf ("-- Time for %4li bits (discriminant %"PRIicl
         ", invariant %c, parameters %s): %.1f\n",
         mpz_sizeinbase (p, 2), d, param->invariant, param->str,
         cm_timer_get (clock));
      if (debug) {
         printf ("   CM:    %5.1f\n", cm_timer_get (stat->timer [1]));
         printf ("   roots: %5.1f\n", cm_timer_get (stat->timer [2]));
         printf ("   point: %5.1f\n", cm_timer_get (stat->timer [3]));
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
   const char* modpoldir, bool tower, bool verbose, bool debug,
   cm_stat_ptr stat)
   /* Given the result of the ECPP down-run in cert1, an array of
      length depth as computed by cm_ecpp1, execute the second step of
      the ECPP algorithm and compute the certificate proper in cert2,
      which needs to be pre-allocated as an array of length depth
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
      the case that verbose is set as well.
      The stat parameter is used to pass timing information to the calling
      function. */

{
   int i;
#ifdef WITH_MPI
   cm_stat_t stat_worker;
   MPI_Status status;
   int sent, received, rank, job;
#endif

   cm_stat_init (stat);
   cm_timer_start (stat->timer [0]);

#ifndef WITH_MPI
   for (i = 0; i < depth; i++)
      cm_ecpp_one_step2 (cert2 [i], cert1 [i], modpoldir, tower,
         verbose, debug, stat);
#else
   sent = 0;
   received = 0;
   while (received < depth) {
      if (sent < depth && (rank = cm_mpi_queue_pop ()) != -1) {
         cm_mpi_submit_ecpp_one_step2 (rank, sent, cert1 [sent], modpoldir,
            tower);
         sent++;
      }
      else {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         cm_mpi_get_ecpp_one_step2 (cert2 [job], rank, stat_worker);
         for (i = 1; i <= 3; i++)
            stat->timer [i]->elapsed += stat_worker->timer [i]->elapsed;
         cm_mpi_queue_push (rank);
         received++;
         if (verbose && debug)
               printf ("-- Timings after job %3i: CM %5.1f, roots %5.1f, "
                  "point %5.1f\n", job,
                  cm_timer_get (stat->timer [1]),
                  cm_timer_get (stat->timer [2]),
                  cm_timer_get (stat->timer [3]));
      }
    }
#endif
   cm_timer_stop (stat->timer [0]);

   if (verbose)
      printf ("--- Time for second ECPP step: %.1f (%.1f)\n",
         cm_timer_get (stat->timer [1]) + cm_timer_get (stat->timer [2])
         + cm_timer_get (stat->timer [3]),
         cm_timer_wc_get (stat->timer [0]));
}

/*****************************************************************************/

bool cm_ecpp (mpz_srcptr N, const char* modpoldir, bool tower,
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
   cm_timer_t clock;
   cm_stat_t stat1, stat2;
   double t;

   cert1 = cm_ecpp1 (&depth, N, verbose, debug, stat1);

   cert2 = (mpz_t **) malloc (depth * sizeof (mpz_t *));
   for (i = 0; i < depth; i++) {
      cert2 [i] = (mpz_t *) malloc (6 * sizeof (mpz_t));
      for (j = 0; j < 6; j++)
         mpz_init (cert2 [i][j]);
   }
   cm_ecpp2 (cert2, cert1, depth, modpoldir, tower, verbose, debug, stat2);
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

   if (verbose) {
      for (i = 0, t = 0.0; i < 7; i++)
         t += cm_timer_get (stat1->timer [i]);
      for (i = 1; i <= 3; i++)
         t += cm_timer_get (stat2->timer [i]);
      printf ("--- Total time for ECPP:       %.1f (%.1f)\n", t,
         cm_timer_wc_get (stat1->timer [7])
         + cm_timer_wc_get (stat2->timer [0]));
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
