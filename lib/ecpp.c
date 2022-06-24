/*

ecpp.c - code for computing ECPP certificates

Copyright (C) 2021, 2022 Andreas Enge

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

static void compute_h (unsigned int *h, uint_cl_t Dmax, cm_stat_t stat);
static void compute_qstar (long int *qstar, mpz_srcptr p, long int *q,
   int no);
static int_cl_t** compute_signed_discriminants (int *no_d, long int *qstar,
   int no_qstar, uint_cl_t Dmax, int sign);
static int_cl_t** compute_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, unsigned int *h);
static int disc_cmp (const void* d1, const void* d2);
static int_cl_t* compute_sorted_discriminants (int *no_d, long int *qstar,
   int no_qstar_old, int no_qstar_new, unsigned int max_factors,
   uint_cl_t Dmax, uint_cl_t hmaxprime, unsigned int *h);
static double expected_no_curves (long int *qstar, int no_qstar_old,
   int no_qstar_new, unsigned int max_factors, uint_cl_t Dmax,
   uint_cl_t hmaxprime, unsigned int *h);
static void compute_qroot (mpz_t *qroot, long int *qstar, int no_qstar,
#ifndef WITH_MPI
   mpz_srcptr p,
#endif
   cm_stat_t stat);
static void sqrt_d (mpz_ptr Droot, int_cl_t d, mpz_srcptr N,
   long int *qstar, int no_qstar, mpz_t *qroot);
static int curve_cardinalities (mpz_t *n, mpz_srcptr N,
   int_cl_t d, long int *qstar, int no_qstar, mpz_t *qroot);
static int card_cmp (const void* c1, const void* c2);
static void trial_div (mpz_t *l, mpz_t *n, int no_n,
#ifndef WITH_MPI
   mpz_srcptr primorialB,
#else
   unsigned long int B,
#endif
   cm_stat_t stat);
static int_cl_t contains_ecpp_discriminant (mpz_ptr n, mpz_ptr l,
   mpz_srcptr N, mpz_t *card, mpz_t *l_list, int_cl_t *d, int no_card,
   const unsigned int delta, cm_stat_t stat);
static int_cl_t find_ecpp_discriminant (mpz_ptr n, mpz_ptr l, mpz_srcptr N,
   uint_cl_t Dmax, uint_cl_t hmaxprime, unsigned int *h,
   const unsigned int delta,
#ifndef WITH_MPI
   mpz_srcptr primorialB,
#endif
   unsigned long int B, bool debug, cm_stat_t stat);
static void ecpp_param_init (cm_param_ptr param, uint_cl_t d);
static mpz_t** ecpp1 (int *depth, mpz_srcptr p, char *filename,
   bool verbose, bool debug, cm_stat_ptr stat);
static void ecpp2 (mpz_t **cert2, mpz_t **cert1, int depth,
   char *filename, const char* modpoldir, bool verbose,
   bool debug, cm_stat_ptr stat);

/*****************************************************************************/

void cm_ecpp_compute_h_chunk (unsigned int *h, uint_cl_t Dmin, uint_cl_t Dmax)
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

static void compute_h (unsigned int *h, uint_cl_t Dmax, cm_stat_t stat)
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
   double t, t_worker;
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
         cm_mpi_get_h_chunk (h + job * size / 2, rank, &t_worker);
         t += t_worker;
         cm_mpi_queue_push (rank);
         received++;
      }
   }
   cm_timer_stop (stat->timer [5]);
   stat->timer [5]->elapsed = t;
#endif
}

/*****************************************************************************/

static void compute_qstar (long int *qstar, mpz_srcptr p, long int *q,
   int no)
   /* Compute and return via qstar an array of no "signed primes" suitable
      for dividing the discriminant to prove the primality of p. The
      entries in qstar are ordered by increasing absolute value (and with
      -8 coming before 8 to break the tie). On first call, *q should be 0;
      it is then replaced by the last found signed prime, and should be
      given for later calls to continue with the next signed primes.
      qstar needs to be initialised to the correct size. */
{
   int i;

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
   bool small;
   int k;

   /* The algorithm assumes that at least the primes in qstar are of
      suitable size; otherwise, forget the largest ones. */
   while (no_qstar > 0
          && (  (qstar [no_qstar-1] > 0
                 && (uint_cl_t) (qstar [no_qstar-1]) > Dmax)
             || (qstar [no_qstar-1] < 0
                 && (uint_cl_t) (-qstar [no_qstar-1]) > Dmax)))
      no_qstar--;

   /* Compute the maximal number of factors that can fit under Dmax. */
   no_factors = 0;
   mpz_init_set_ui (tmp, 1);
   small = true;
   while (no_factors < no_qstar && small) {
      if (qstar [no_factors] > 0)
         mpz_mul_ui (tmp, tmp, (unsigned long int) (qstar [no_factors]));
      else
         mpz_mul_ui (tmp, tmp, (unsigned long int) (-qstar [no_factors]));
      if (mpz_cmp_ui (tmp, Dmax) <= 0)
         no_factors++;
      else
         small = false;
   }

   /* Compute an upper bound on the possible number of discriminants. */
   mpz_init (bin);
   mpz_set_ui (tmp, 0);
   for (k = 0; k <= no_factors; k++) {
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
   if (sign == 1 && Dmax > 0) {
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
               need to go one level up. Notice that if -8 and 8 are elements
               of qstar, the list is not strictly increasing: So when
               Dmax==8, we have to continue even when |D|==Dmax. */
            if (Dqmax == no_qstar - 1
                  || (D > 0 && (uint_cl_t)   D  > Dmax)
                  || (D < 0 && (uint_cl_t) (-D) > Dmax)) {
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
   uint_cl_t Dmax, uint_cl_t hmaxprime, unsigned int *h)
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
   int i, j, k;
   long int exclude[50];
   int no_exclude = 0;
   FILE *exclude_file;

   exclude_file = fopen("exclude.txt", "r");
   if (exclude_file != NULL) {
      while (fscanf(exclude_file, "%ld", &exclude[no_exclude]) == 1) {
         no_exclude++;
         if (no_exclude == 50) break;
      }
      fclose(exclude_file);
   }

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
         for (k = 0; k < no_exclude; k++)
            if (D == (int_cl_t)exclude[k]) D = 16; 
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
   uint_cl_t Dmax, uint_cl_t hmaxprime, unsigned int *h)
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
   if (*no_d > 0)
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

static double expected_no_curves (long int *qstar, int no_qstar_old,
   int no_qstar_new, unsigned int max_factors, uint_cl_t Dmax,
   uint_cl_t hmaxprime, unsigned int *h)
   /* The function takes the same parameters as compute_discriminants (and
      most of them are just passed through). It returns the expected
      number of curve cardinalities obtained from the list of discriminants
      computed by compute_discriminants. */
{
   int_cl_t **dlist;
   int no_d;
   double exp_card;
      /* the expected number of curve cardinalities; the sum over
         #twists * (g/h) */
   int i;

   dlist = compute_discriminants (&no_d, qstar, no_qstar_old, no_qstar_new,
      max_factors, Dmax, hmaxprime, h);

   exp_card = 0;
   for (i = 0; i < no_d; i++) {
      if (dlist [i][0] == -3)
         exp_card += 6.0;
      else if (dlist [i][0] == -4)
         exp_card += 4.0;
      else
         exp_card += 2.0 / dlist [i][2];
      free (dlist [i]);
   }
   free (dlist);

   return exp_card;
}

/*****************************************************************************/

static void compute_qroot (mpz_t *qroot, long int *qstar, int no_qstar,
#ifndef WITH_MPI
   mpz_srcptr p,
#endif
   cm_stat_t stat)
   /* Compute in qroot the square roots of the no_qstar elements of qstar
      (assumed to be squares) modulo p. qroot needs to be initialised to
      the correct size, and its entries need to be initialised as well. */
{
#ifndef WITH_MPI
   unsigned int e;
   mpz_t r, z;
   int i;
#else
   MPI_Status status;
   int sent, received, rank, job;
   double t, t_worker;
#endif

#ifndef WITH_MPI
   cm_timer_continue (stat->timer [0]);
   mpz_init (r);
   mpz_init (z);
   e = cm_nt_mpz_tonelli_generator (r, z, p);
   for (i = 0; i < no_qstar; i++)
      cm_nt_mpz_tonelli_si_with_generator (qroot [i], qstar [i],
         p, e, r, z);
   mpz_clear (r);
   mpz_clear (z);
   cm_timer_stop (stat->timer [0]);
#else
   sent = 0;
   received = 0;
   t = cm_timer_get (stat->timer [0]);
      /* Memorise CPU time to avoid confusion with server CPU time. */
   cm_timer_continue (stat->timer [0]);
   while (received < no_qstar) {
      if (sent < no_qstar && (rank = cm_mpi_queue_pop ()) != -1) {
         cm_mpi_submit_tonelli (rank, sent, qstar [sent]);
         sent++;
      }
      else {
         /* There is either nothing to send, or all workers are busy.
            Wait for a result. */
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         cm_mpi_get_tonelli (qroot [job], rank, &t_worker);
         t += t_worker;
         cm_mpi_queue_push (rank);
         received++;
      }
   }
   cm_timer_stop (stat->timer [0]);
   stat->timer [0]->elapsed = t;
      /* Restore CPU time containing all the workers, while keeping
         the wallclock time as measured by the server. */
#endif
}

/*****************************************************************************/

static void sqrt_d (mpz_ptr Droot, int_cl_t d, mpz_srcptr N,
   long int *qstar, int no_qstar, mpz_t *qroot)
   /* Given a discriminant d which factors over the array of no_qstar
      signed primes qstar, and given the roots of qstar modulo N in qroot,
      compute and return in Droot the root of the d.
      By trial dividing d over qstar, it is enough to multiply the
      corresponding roots together; however, the "even primes" need special
      care, in particular for the sign of 8. */
{
   int i;

   mpz_set_ui (Droot, 1);
   for (i = 0; i < no_qstar; i++)
      if (qstar [i] % 2 != 0 && d % qstar [i] == 0) {
         mpz_mul (Droot, Droot, qroot [i]);
         mpz_mod (Droot, Droot, N);
         d /= qstar [i];
      }
   if (d != 1) {
      for (i = 0; d != qstar [i]; i++);
      mpz_mul (Droot, Droot, qroot [i]);
      mpz_mod (Droot, Droot, N);
   }
}

/*****************************************************************************/

static int curve_cardinalities (mpz_t *n, mpz_srcptr N, int_cl_t d,
   long int *qstar, int no_qstar, mpz_t *qroot)
   /* Given N prime, a list of no_qstar signed primes in qstar, their
      roots modulo N in qroot and a discriminant d that factors over qstar,
      the function computes the array of 0 (in the case that N is not a
      norm in Q(\sqrt d), which happens with probability 1 - 1 / (h/g)),
      2, 4 or 6 (depending on the number of twists) possible cardinalities
      of elliptic curves modulo N with CM by d. The cardinalities are
      stored in n, the entries of which need to be initialised, and their
      number is returned. */
{
   int res;
   mpz_t root, t, v;
   mpz_ptr V;
   int twists;
   bool cornacchia;

   mpz_init (root);
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

   sqrt_d (root, d, N, qstar, no_qstar, qroot);
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
   mpz_clear (root);
   mpz_clear (t);

   return res;
}

/*****************************************************************************/

extern mpz_t* cm_ecpp_compute_cardinalities (int *no_card,
   int_cl_t **card_d, int_cl_t *d, int no_d, mpz_srcptr N,
   long int *qstar, int no_qstar, mpz_t *qroot)
   /* Given a prime N, a list of no_d fastECPP discriminants in d factoring
      over the array qstar of no_qstar signed primes and the array qroot
      of the square roots of qstar modulo N, compute and return the
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
      twists [i] = curve_cardinalities (card [i], N, d [i],
         qstar, no_qstar, qroot);

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

void trial_div (mpz_t *l, mpz_t *n, int no_n,
#ifndef WITH_MPI
   mpz_srcptr primorialB,
#else
   unsigned long int B,
#endif
   cm_stat_t stat)
   /* primorialB is supposed to be the product of the primes up to the
      smoothness bound B.
      The function removes all occurrences of the primes up to B from the
      no_n entries in n and stores the results in l, which needs to have
      the correct size and all entries of which need to be initialised. */
{
   mpz_t *gcd;
   int i;
#ifdef WITH_MPI
   double t;
#endif

   cm_timer_continue (stat->timer [2]);
   gcd = (mpz_t *) malloc (no_n * sizeof (mpz_t));
   for (i = 0; i < no_n; i++)
      mpz_init (gcd [i]);

   /* Compute in gcd [i] the gcd of n [i] and primorialB. */
#ifndef WITH_MPI
   cm_nt_mpz_tree_gcd (gcd, primorialB, n, no_n);
#else
   cm_mpi_submit_tree_gcd (n, no_n);
   cm_mpi_get_tree_gcd (gcd, no_n, B, &t);
#endif
   cm_timer_stop (stat->timer [2]);
   stat->counter [2] += no_n;
#ifdef WITH_MPI
   stat->timer [2]->elapsed += t;
#endif

   cm_timer_continue (stat->timer [2]);
   /* Remove the gcd from n [i] and recompute it until all primes
      are removed. */
   for (i = 0; i < no_n; i++) {
      mpz_set (l [i], n [i]);
      while (mpz_cmp_ui (gcd [i], 1ul) != 0) {
         mpz_divexact (l [i], l [i], gcd [i]);
         mpz_gcd (gcd [i], l [i], gcd [i]);
      }
   }

   for (i = 0; i < no_n; i++)
      mpz_clear (gcd [i]);
   free (gcd);
   cm_timer_stop (stat->timer [2]);
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
   int no;
   mpz_t **c;
   int i;
#ifndef WITH_MPI
   int index;
#else
   int *index;
   int batch, max_i;
   int j;
   MPI_Status status;
   int size, rank, job;
   double t, t_worker;
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
   if (no > 0)
      qsort (c, no, sizeof (mpz_t *), card_cmp);

#ifndef WITH_MPI
   /* Find the first prime non-smooth part, which is also the
      smallest one. */
   cm_timer_continue (stat->timer [3]);
   for (i = 0; res == 0 && i < no; i++) {
      if (cm_nt_is_prime (c [i][0])) {
         index = mpz_get_ui (c [i][1]);
         res = d [index];
         mpz_set (n, card [index]);
         mpz_set (l, l_list [index]);
      }
      stat->counter [3]++;
   }
   cm_timer_stop (stat->timer [3]);
#else
   /* Let each worker do one primality test at a time; stop as soon as one
      of the workers finds a prime, then keep the first one in case several
      workers identify a prime in one batch. */
   t = cm_timer_get (stat->timer [3]);
   cm_timer_continue (stat->timer [3]);
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   batch = size - 1;
   index = (int *) malloc (batch * sizeof (int));
   for (j = 0; j < batch; j++)
      index [j] = -1;
   max_i = (no + batch - 1) / batch;
   for (i = 0; res == 0 && i < max_i; i++) {
      if (i == max_i - 1)
         batch = no - i * batch;
      for (j = 0; j < batch; j++)
         cm_mpi_submit_is_prime (j + 1, j, c [j + i * batch][0]);
      for (j = 0; j < batch; j++) {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         if (cm_mpi_get_is_prime (rank, &t_worker))
            index [job] = mpz_get_ui (c [job + i * batch][1]);
         t += t_worker;
      }
      for (j = 0; res == 0 && j < batch; j++)
         if (index [j] != -1) {
            res = d [index [j]];
            mpz_set (n, card [index [j]]);
            mpz_set (l, l_list [index [j]]);
         }
      stat->counter [3] += batch;
   }
   cm_timer_stop (stat->timer [3]);
   stat->timer [3]->elapsed = t;
   free (index);
#endif

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
   uint_cl_t Dmax, uint_cl_t hmaxprime, unsigned int *h,
   const unsigned int delta,
#ifndef WITH_MPI
   mpz_srcptr primorialB,
#endif
   unsigned long int B, bool debug, cm_stat_t stat)
   /* Given a (probable) prime N>=787, return a suitable CM discriminant
      and return the cardinality of an associated elliptic curve in n and
      its largest prime factor in l.
      Dmax, hmaxprime and h are passed through to compute_discriminants.
      delta >= 1 is passed through as the minimum number of bits to be
      gained in this step.
      primorialB is passed through to trial division.
      The smoothness bound B is needed to compute the success probability
      of finding a prime. */
{
   const int max_factors = 7;
   int no_qstar_old, no_qstar_new, no_qstar, no_qstar_delta;
   long int *qstar;
   long int q;
   mpz_t *root, *card, *l_list;
   int_cl_t d;
   int_cl_t *dlist, *d_card;
   int no_d, no_card;
   const double prob_prime = 1.7811
      * log2 (B) / mpz_sizeinbase (N, 2);
      /* According to [FrKlMoWi04], the probability that a number N is
         a B-smooth part times a prime is exp (gamma) * log B / log N. */
   double exp_prime, min_prime;
   int round, i;
#ifdef WITH_MPI
   MPI_Status status;
   int size, rank, job;
   double t, t_worker;
   int batch, max_i;
   int_cl_t **d_card_batch;
   mpz_t **card_batch;
   int *no_card_batch;
   int j;
#endif

   d = 0;
   no_qstar = 0;
   q = 0;
   qstar = (long int *) malloc (0);
   root = (mpz_t *) malloc (0);

#ifdef WITH_MPI
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   cm_mpi_broadcast_N (N);
   no_qstar_delta = size - 1;
#else
   no_qstar_delta = 1;
#endif
   min_prime = 3.0;
      /* With this expected number of suitable cardinalities, there
         should be a good chance of finding at least one; more lead
         to some choice for gaining more bits in one step. */
   round = 0;
   while (d == 0) {
      round++;
      if (debug)
         printf ("  Round %i\n  no_qstar_delta: %i\n",
            round, no_qstar_delta);
      /* Extend the prime list in small pieces until the expectation of
         finding a curve cardinality that is a smooth part times a prime
         becomes sufficiently large. */
      cm_timer_continue (stat->timer [4]);
      no_qstar_old = no_qstar;
      exp_prime = 0.0;
      while (exp_prime < min_prime) {
         no_qstar += no_qstar_delta;
         qstar = (long int *) realloc (qstar, no_qstar * sizeof (long int));
         compute_qstar (qstar + (no_qstar - no_qstar_delta), N, &q,
            no_qstar_delta);
         exp_prime += prob_prime *
            expected_no_curves (qstar, no_qstar - no_qstar_delta,
               no_qstar_delta, max_factors, Dmax, hmaxprime, h);
         if (debug)
            printf ("    no_qstar: %i, exp_prime: %.1f\n",
               no_qstar, exp_prime);
      }

      /* Now compute the list of discriminants containing one of the
         new primes. */
      no_qstar_new = no_qstar - no_qstar_old;
      dlist = compute_sorted_discriminants (&no_d, qstar, no_qstar_old,
         no_qstar_new, max_factors, Dmax, hmaxprime, h);
      cm_timer_stop (stat->timer [4]);
      if (debug)
         printf ("    no_d: %i\n", no_d);

      /* Compute (and broadcast in the case of MPI) the square roots of
         the new primes. */
      root = (mpz_t *) realloc (root, no_qstar * sizeof (mpz_t));
      for (i = no_qstar_old; i < no_qstar; i++)
         mpz_init (root [i]);
      compute_qroot (root + no_qstar_old, qstar + no_qstar_old,
         no_qstar_new,
#ifndef WITH_MPI
         N,
#endif
         stat);
#ifdef WITH_MPI
      cm_timer_continue (stat->timer [0]);
      cm_mpi_broadcast_sqrt (no_qstar_new, qstar + no_qstar_old,
         root + no_qstar_old);
      cm_timer_stop (stat->timer [0]);
#endif
      stat->counter [0] += no_qstar_new;

      /* Compute the cardinalities of the corresponding elliptic curves. */
#ifndef WITH_MPI
      cm_timer_continue (stat->timer [1]);
      card = cm_ecpp_compute_cardinalities (&no_card, &d_card, dlist, no_d,
         N, qstar, no_qstar, root);
      cm_timer_stop (stat->timer [1]);
      stat->counter [1] += no_d;
#else
      cm_timer_continue (stat->timer [1]);
      batch = (no_d + size - 2) / (size - 1);
      max_i = (no_d + batch - 1) / batch;
      no_card_batch = (int *) malloc (max_i * sizeof (int));
      d_card_batch = (int_cl_t **) malloc (max_i * sizeof (int_cl_t *));
      card_batch = (mpz_t **) malloc (max_i * sizeof (mpz_t *));
      for (i = 0; i < max_i; i++) {
         rank = cm_mpi_queue_pop ();
         cm_mpi_submit_curve_cardinalities (rank, i, dlist + i * batch,
            (i == max_i - 1 ? no_d - (max_i - 1) * batch : batch));
      }
      cm_timer_stop (stat->timer [1]);
      t = cm_timer_get (stat->timer [1]);
      cm_timer_continue (stat->timer [1]);
      for (i = 0; i < max_i; i++) {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         card_batch [i] = cm_mpi_get_curve_cardinalities (no_card_batch + i,
            d_card_batch + i, rank, &t_worker);
         t += t_worker;
         cm_mpi_queue_push (rank);
      }
      for (no_card = 0, i = 0; i < max_i; i++)
         no_card += no_card_batch [i];
      d_card = (int_cl_t *) malloc (no_card * sizeof (int_cl_t));
      card = (mpz_t *) malloc (no_card * sizeof (mpz_t));
      for (no_card = 0, i = 0; i < max_i; i++) {
         for (j = 0; j < no_card_batch [i]; j++) {
            d_card [no_card + j] = d_card_batch [i][j];
            card [no_card + j][0] = card_batch [i][j][0];
         }
         no_card += no_card_batch [i];
         free (d_card_batch [i]);
         free (card_batch [i]);
      }
      free (no_card_batch);
      free (d_card_batch);
      free (card_batch);
      cm_timer_stop (stat->timer [1]);
      stat->timer [1]->elapsed = t;
      stat->counter [1] += no_d;
#endif

      if (no_card > 0) {
         /* Remove smooth parts of cardinalities. */
         l_list = (mpz_t *) malloc (no_card * sizeof (mpz_t));
         for (i = 0; i < no_card; i++)
            mpz_init (l_list [i]);
         trial_div (l_list, card, no_card,
#ifndef WITH_MPI
            primorialB,
#else
            B,
#endif
            stat);

         d = contains_ecpp_discriminant (n, l, N, card, l_list, d_card,
               no_card, delta, stat);

         for (i = 0; i < no_card; i++)
            mpz_clear (l_list [i]);
         free (l_list);
      }

      free (dlist);
      for (i = 0; i < no_card; i++)
         mpz_clear (card [i]);
      free (card);
      free (d_card);

      /* If no discriminants are found in the first round, chances are that
         the class numbers have become quite high, and that the expected
         number of curve cardinalities per discriminant quite low; so we
         should lower our expectations. */
      if (min_prime == 3.0)
         min_prime = 1.0;
      else if (min_prime == 1.0)
         min_prime = 0.5;
      else
         min_prime = 0.25;
   }

#ifdef WITH_MPI
      cm_mpi_clear_N ();
#endif
   for (i = 0; i < no_qstar; i++)
      mpz_clear (root [i]);
   free (root);
   free (qstar);

   return d;
}

/*****************************************************************************/

static mpz_t** ecpp1 (int *depth, mpz_srcptr p, char *filename,
   bool verbose, bool debug, cm_stat_ptr stat)
   /* Compute the first step of the ECPP certificate; this is the downrun
      part with the parameters of the elliptic curves.
      The return value is a newly allocated array of depth entries, each
      of which is an array of length 4, containing in this order
      - p_i, a prime to be certified;
      - d_i, the discriminant;
      - n_i, the cardinality of the elliptic curve;
      - l_i, the prime order dividing this cardinality.
      The downrun stops as soon as the prime is less than 2^64.
      If filename is different from NULL, then the function tries to read
      a partial certificate from the file and continues writing to it.
      The stat parameter makes it possible to pass timing information
      to the calling function. */

{
   const size_t L = mpz_sizeinbase (p, 2);
#ifndef WITH_MPI
   const unsigned long int B = (L >= 26008 ? 1ul<<31 :
                                             (L>>4)*(L>>4)*(L>>5));
   const unsigned int delta = (unsigned int) (log2 (B) / 2) + 1;
      /* According to [FrKlMoWi04] the average factor removed by trial
         division up to B, assuming that what remains is prime, is B;
         we impose half of this number of bits as the minimal gain. */
#else
   unsigned long int B;
      /* Computed below depending on size. */
   const unsigned int delta = 2;
      /* Since in the parallel version we consider many potential curve
         orders at once and order them by gain, having a small value of
         delta does not make much difference most of the time. For "bad
         primes", it enables us to side-step the issue instead of adding
         more and more qstar, which becomes increasingly difficult. */
#endif
   const uint_cl_t Dmax = ((L * L) >> 4) << 2;
   const uint_cl_t hmaxprime = CM_MAX (29, L>>10);
   mpz_t N;
   mpz_t** c;
   unsigned int *h;
   int_cl_t d;
   FILE *f;
   cm_timer_t clock;
   double t, t_old;
   int i;
#ifndef WITH_MPI
   mpz_t primorialB;
#else
   int size, rank, job;
   MPI_Status status;
   double t_worker;
#endif

#ifdef WITH_MPI
   /* Cap B at (size - 1) * 2^29, so that each worker handles a product
      of value at most about exp (2^29), or 100MB. */
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   B = CM_MIN (((unsigned long int) size - 1) << 29, (L>>4)*(L>>4)*(L>>5));
#endif

   cm_stat_init (stat);

   /* If filename is given, try to read a partial step 1 certificate from
      the file. */
   mpz_init (N);
   *depth = 0;
   if (filename != NULL)
      if (cm_file_open_read_write (&f, filename))
         c = cm_file_read_ecpp_cert1 (depth, p, f, debug, stat);
      else {
         cm_file_open_write (&f, filename);
         c = NULL;
      }
   else
      c = NULL;
   if (*depth > 0)
      mpz_set (N, c [*depth - 1][3]);
   else
      mpz_set (N, p);

   cm_timer_continue (stat->timer [7]);

   if (mpz_sizeinbase (N, 2) > 64) {
      /* Precompute class numbers. */
      h = (unsigned int *) malloc ((Dmax / 2) * sizeof (unsigned int));
      compute_h (h, Dmax, stat);
      if (verbose)
         printf ("Time for class numbers up to Dmax=%"PRIucl
               ": %.0f (%.0f)\n", Dmax,
               cm_timer_get (stat->timer [5]),
               cm_timer_wc_get (stat->timer [5]));

      /* Precompute primorial for trial division. */
#ifndef WITH_MPI
      mpz_init (primorialB);
      cm_timer_start (stat->timer [6]);
      mpz_primorial_ui (primorialB, B);
      cm_timer_stop (stat->timer [6]);
#else
      t = 0;
      cm_timer_start (stat->timer [6]);
      cm_mpi_submit_primorial (B);
      for (i = 1; i < size; i++) {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status);
         rank = status.MPI_SOURCE;
         cm_mpi_get_primorial (rank, &t_worker);
         t += t_worker;
      }
      cm_timer_stop (stat->timer [6]);
      stat->timer [6]->elapsed = t;
#endif
      if (verbose) {
         printf ("Time for primorial of B=%lu: %.0f (%.0f)\n", B,
               cm_timer_get (stat->timer [6]), cm_timer_wc_get (stat->timer [6]));
         printf ("hmaxprime: %"PRIucl"\n", hmaxprime);
      }

      t_old = 0;
      while (mpz_sizeinbase (N, 2) > 64) {
         c = (mpz_t**) realloc (c, (*depth + 1) * sizeof (mpz_t *));
         c [*depth] = (mpz_t *) malloc (4 * sizeof (mpz_t));
         for (i = 0; i < 4; i++)
            mpz_init (c [*depth][i]);
         mpz_set (c [*depth][0], N);
         cm_timer_start (clock);
         if (verbose)
            printf ("Size [%4i]: %6li bits\n", *depth,
                  mpz_sizeinbase (N, 2));
         d = find_ecpp_discriminant (c [*depth][2], c [*depth][3], N, Dmax,
               hmaxprime, h, delta,
#ifndef WITH_MPI
               primorialB,
#endif
               B, debug, stat);
         cm_timer_stop (clock);
         if (verbose) {
            for (i = 0, t = 0.0; i < 5; i++)
               t += cm_timer_get (stat->timer [i]);
            printf ("  Time for discriminant %11"PRIicl": %.0f (%.0f)\n",
                  d, t - t_old, cm_timer_wc_get (clock));
            t_old = t;
            if (verbose) {
               printf ("  largest prime of d: %"PRIucl"\n",
                     cm_nt_largest_factor (-d));
               printf ("  largest prime of h: %"PRIucl"\n",
                     cm_nt_largest_factor (h [(-d) / 2 - 1]));
               printf ("  discriminants:          %11.0f (%7.0f)\n",
                     cm_timer_get (stat->timer [4]),
                     cm_timer_wc_get (stat->timer [4]));
               printf ("%13lu qroot:      %11.0f (%7.0f)\n", stat->counter [0],
                     cm_timer_get (stat->timer [0]),
                     cm_timer_wc_get (stat->timer [0]));
               printf ("%13lu Cornacchia: %11.0f (%7.0f)\n", stat->counter [1],
                     cm_timer_get (stat->timer [1]),
                     cm_timer_wc_get (stat->timer [1]));
               printf ("%13lu trial div:  %11.0f (%7.0f)\n", stat->counter [2],
                     cm_timer_get (stat->timer [2]),
                     cm_timer_wc_get (stat->timer [2]));
               printf ("%13lu is_prime:   %11.0f (%7.0f)\n", stat->counter [3],
                     cm_timer_get (stat->timer [3]),
                     cm_timer_wc_get (stat->timer [3]));
            }
         }
         mpz_set_si (c [*depth][1], d);
         if (filename != NULL) {
            cm_timer_stop (stat->timer [7]);
            cm_write_ecpp_cert1_line (f, c [*depth], stat);
            cm_timer_continue (stat->timer [7]);
         }
         mpz_set (N, c [*depth][3]);
         (*depth)++;
      }
      free (h);
#ifndef WITH_MPI
      mpz_clear (primorialB);
#endif
   }

   mpz_clear (N);

   cm_timer_stop (stat->timer [7]);
   for (i = 0, t = 0.0; i < 7; i++)
      t+= cm_timer_get (stat->timer [i]);
   if (verbose)
      printf ("Time for first ECPP step, depth %i: %.0f (%.0f)\n",
         *depth, t, cm_timer_wc_get (stat->timer [7]));

   if (filename != NULL)
      cm_file_close (f);

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

void cm_ecpp_one_step2 (mpz_t *cert2, mpz_t *cert1, int i,
   const char* modpoldir, bool verbose, bool debug, cm_stat_t stat)
   /* The function takes the same parameters as ecpp2, except that cert1
      contains only one entry of the first ECPP step, and that only one
      step of the certificate is output in cert2, which needs be to a
      pre-allocated and initialised array with 6 entries. The integer i
      is the number of the current step; it is used only for printing
      progress information.
      This function has been separated to help with parallelisation. */
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

   /* Compute the class field tower. */
   cm_class_init (c, param, false);
   cm_class_compute (c, param, false, true, false);
   cm_timer_stop (stat->timer [1]);

   cm_curve_and_point_stat (a, b, x, y, param, c, p, l, co,
      modpoldir, false, false, stat);
   cm_class_clear (c);
   cm_timer_stop (clock);

   if (verbose) {
      printf ("Time for %4i (discriminant %11"PRIicl
         ", invariant %c, parameters %13s): %6.0f\n",
         i, d, param->invariant, param->str, cm_timer_get (clock));
      if (debug) {
         printf ("  CM:    %5.0f\n", cm_timer_get (stat->timer [1]));
         printf ("  roots: %5.0f\n", cm_timer_get (stat->timer [2]));
         printf ("  point: %5.0f\n", cm_timer_get (stat->timer [3]));
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

static void ecpp2 (mpz_t **cert2, mpz_t **cert1, int depth, char *filename,
   const char* modpoldir, bool verbose, bool debug, cm_stat_ptr stat)
   /* Given the result of the ECPP down-run in cert1, an array of
      length depth as computed by ecpp1, execute the second step of
      the ECPP algorithm and compute the certificate proper in cert2,
      which needs to be pre-allocated as an array of length depth
      with each entry an array of length 6 to contain p, t, co, a, x and y:
      the prime p to be certified, the trace t of the elliptic curve, the
      co-factor such that p+1-t = co*l with l the next prime, the curve
      parameter a, and the coordinates (x,y) of a point of prime order l
      on the curve; the curve parameter b is implicit.
      If different from NULL, filename gives a file from which to read the
      beginning of the certificate, and to which to write the computed end.
      modpoldir gives the directory where modular polynomials are stored;
      it is passed through to the function computing a curve from a root
      of the class polynomial.
      verbose indicates whether intermediate computations output progress
      information.
      debug indicates whether additional developer information (mainly
      timings and counters for tuning) is output; this is done only in
      the case that verbose is set as well.
      The stat parameter is used to pass timing information to the calling
      function. */

{
   int done;
   int i;
   FILE *f;
#ifdef WITH_MPI
   cm_stat_t stat_worker;
   MPI_Status status;
   int sent, received, rank, job;
#endif

   cm_stat_init (stat);

   if (filename != NULL)
      if (cm_file_open_read_write (&f, filename))
         cm_file_read_ecpp_cert2 (cert2, cert1 [0][0], f, debug,
            stat);
      else {
         cm_file_open_write (&f, filename);
      }

   cm_timer_continue (stat->timer [0]);
#ifndef WITH_MPI
   for (i = 0; i < depth; i++) {
      if (mpz_sgn (cert2 [i][0]) == 0) {
         cm_ecpp_one_step2 (cert2 [i], cert1 [i], i, modpoldir,
            verbose, debug, stat);
         if (filename != NULL) {
            cm_timer_stop (stat->timer [0]);
            cm_write_ecpp_cert2_line (f, i, cert2 [i], stat);
            cm_timer_continue (stat->timer [0]);
         }
      }
   }
#else
   sent = 0;
   received = 0;
   for (i = 0; i < depth; i++) 
      if (mpz_sgn (cert2 [i][0]) != 0) received++;
   
   while (received < depth) {
      while (sent < depth && mpz_sgn (cert2 [sent][0]) != 0) sent++;
      if (sent < depth && (rank = cm_mpi_queue_pop ()) != -1) {
         cm_mpi_submit_ecpp_one_step2 (rank, sent, cert1 [sent], modpoldir);
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
         if (debug)
            printf ("Timings after job %3i: CM %7.0f, roots %10.0f, "
               "point %7.0f\n", job,
               cm_timer_get (stat->timer [1]),
               cm_timer_get (stat->timer [2]),
               cm_timer_get (stat->timer [3]));
         if (filename != NULL) {
            cm_timer_stop (stat->timer [0]);
            cm_write_ecpp_cert2_line (f, job, cert2 [job], stat);
            cm_timer_continue (stat->timer [0]);
         }
      }
    }
#endif
   cm_timer_stop (stat->timer [0]);

   if (verbose)
      printf ("Time for second ECPP step: %.0f (%.0f)\n",
         cm_timer_get (stat->timer [1]) + cm_timer_get (stat->timer [2])
         + cm_timer_get (stat->timer [3]),
         cm_timer_wc_get (stat->timer [0]));
}

/*****************************************************************************/

bool cm_ecpp (mpz_srcptr N, const char* modpoldir, bool print,
   char *filename, bool trust, bool check, bool verbose, bool debug, bool onlys1)
   /* Assuming that N is a (probable) prime, compute an ECPP certificate.
      modpoldir gives the directory where modular polynomials are stored;
      it is passed through to the function computing a curve from a root
      of the class polynomial.
      If filename is different from NULL, the final ECPP certificate is
      output to the file, and the stage 1 and stage 2 certificates are
      read from (partially) and written to temporary files.
      If trust is set to true, then N is trusted to be a probable prime;
      otherwise a quick primality test is run.
      print indicates whether the result is printed to stdout.
      check indicates whether the certificate should be checked.
      If yes, the return value of the function is the result of the check;
      otherwise the return value is true.
      verbose indicates whether intermediate computations output progress
      information.
      debug indicates whether additional developer information (mainly
      timings and counters for tuning) is output; this is done only in
      the case that verbose is set as well. 
      onlys1 indicates to stop at the end of stage 1*/

{
   bool res = true;
   int depth;
   mpz_t **cert1, **cert2;
   char *filename1, *filename2, *filenameprimo;
   int i, j;
   cm_timer_t clock;
   cm_stat_t stat1, stat2;
   double t;
   FILE *f;

   if (!trust) {
      cm_timer_start (clock);
      if (!mpz_probab_prime_p (N, 1)) {
         printf ("***** Error: cm_ecpp called with composite number.\n");
         exit (1);
      }
      else {
         cm_timer_stop (clock);
         if (verbose)
            printf ("Time for primality test: %.0f (%.0f)\n",
            cm_timer_get (clock), cm_timer_wc_get (clock));
      }
   }

   if (filename != NULL) {
      i = strlen (filename) + 7;
      filename1 = (char *) malloc (i * sizeof (char));
      sprintf (filename1, "%s.cert1", filename);
      filename2 = (char *) malloc (i * sizeof (char));
      sprintf (filename2, "%s.cert2", filename);
      filenameprimo = (char *) malloc (i * sizeof (char));
      sprintf (filenameprimo, "%s.primo", filename);
   }
   else {
      filename1 = NULL;
      filename2 = NULL;
   }
   cert1 = ecpp1 (&depth, N, filename1, verbose, debug, stat1);

   if (!onlys1) {
      cert2 = (mpz_t **)malloc(depth * sizeof(mpz_t *));
      for (i = 0; i < depth; i++)
      {
         cert2[i] = (mpz_t *)malloc(6 * sizeof(mpz_t));
         for (j = 0; j < 6; j++)
            mpz_init(cert2[i][j]);
      }
      ecpp2(cert2, cert1, depth, filename2, modpoldir, verbose, debug, stat2);

      if (print)
         cm_file_write_ecpp_cert_pari(stdout, cert2, depth);

      if (filename != NULL)
      {
         if (!cm_file_open_write(&f, filename))
            exit(1);
         cm_file_write_ecpp_cert_pari(f, cert2, depth);
         cm_file_close(f);
         if (!cm_file_open_write(&f, filenameprimo))
            exit(1);
         cm_file_write_ecpp_cert_primo(f, cert2, depth);
         cm_file_close(f);
      }

      if (verbose)
      {
         for (i = 0, t = 0.0; i <= 6; i++)
            t += cm_timer_get(stat1->timer[i]);
         for (i = 1; i <= 3; i++)
            t += cm_timer_get(stat2->timer[i]);
         printf("Total time for ECPP: %.0f (%.0f)\n", t,
                cm_timer_wc_get(stat1->timer[7]) + cm_timer_wc_get(stat2->timer[0]));
      }

      if (check)
      {
         cm_timer_start(clock);
         res = cm_pari_ecpp_check(cert2, depth);
         cm_timer_stop(clock);
         if (verbose)
            printf("Time for ECPP check (%s): %.1f\n",
                   (res ? "true" : "false"), cm_timer_get(clock));
      }
   }

   for (i = 0; i < depth; i++) {
      for (j = 0; j < 4; j++)
         mpz_clear (cert1 [i][j]);
      free (cert1 [i]);
      if (!onlys1) {
         for (j = 0; j < 6; j++)
            mpz_clear (cert2 [i][j]);
         free (cert2 [i]);
      }
   }
   free (cert1);
   if (!onlys1) free (cert2);
   if (filename != NULL) {
      free (filename1);
      free (filename2);
   }

   return res;
}

/*****************************************************************************/
/*****************************************************************************/
