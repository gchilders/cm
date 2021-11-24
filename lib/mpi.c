/*

mpicm.c - functions enabling MPI for CM

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

#define MPI_TAG_FINISH                   1
#define MPI_TAG_DATA                     2
#define MPI_TAG_JOB_TONELLI              3
#define MPI_TAG_JOB_ECPP2                4
#define MPI_TAG_JOB_CARD                 5
#define MPI_TAG_JOB_PRIME                6
#define MPI_TAG_JOB_H                    7
#define MPI_TAG_JOB_TRIAL_DIV            8
#define MPI_TAG_JOB_BROADCAST_N          9
#define MPI_TAG_JOB_BROADCAST_PRIMORIAL 10
#define MPI_TAG_JOB_BROADCAST_SQRT      11
#define MPI_TAG_JOB_CLEAR_N             12

typedef char mpi_name_t [MPI_MAX_PROCESSOR_NAME];

static int *worker_queue, *worker_queue_local;
static int worker_queue_size, worker_queue_local_size;

static void mpi_send_mpz (mpz_srcptr z, const int rank);
static void mpi_recv_mpz (mpz_ptr z, const int rank);
static void mpi_bcast_send_mpz (mpz_srcptr z);
static void mpi_bcast_recv_mpz (mpz_ptr z);
static void mpi_worker (void);
static void mpi_server_init (const int size, bool debug);
static void mpi_server_clear (const int size);

/*****************************************************************************/
/*                                                                           */
/* Data structures sending and receiving.                                    */
/*                                                                           */
/*****************************************************************************/

static void mpi_send_mpz (mpz_srcptr z, const int rank)
   /* Send z to rank. */
{
   int size = z->_mp_size;

   MPI_Send (&size, 1, MPI_INT, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD);
   if (size > 0)
      MPI_Send (z->_mp_d,  size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         MPI_COMM_WORLD);
   else if (size < 0)
      MPI_Send (z->_mp_d, -size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         MPI_COMM_WORLD);
}

/*****************************************************************************/

static void mpi_recv_mpz (mpz_ptr z, const int rank)
   /* Get z from rank. */
{
   int size;

   MPI_Recv (&size, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
   if (size > 0) {
      _mpz_realloc (z, size);
      MPI_Recv (z->_mp_d,  size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         MPI_COMM_WORLD, NULL);
   }
   else if (size < 0) {
      _mpz_realloc (z, -size);
      MPI_Recv (z->_mp_d, -size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         MPI_COMM_WORLD, NULL);
   }
   z->_mp_size = size;
}

/*****************************************************************************/

static void mpi_bcast_send_mpz (mpz_srcptr z)
   /* Upon a call by rank 0, send z by broadcast to all others. */
{
   int size = z->_mp_size;

   MPI_Bcast (&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (size > 0)
      MPI_Bcast (z->_mp_d,  size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
   else if (size < 0)
      MPI_Bcast (z->_mp_d, -size, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
}

/*****************************************************************************/

static void mpi_bcast_recv_mpz (mpz_ptr z)
   /* Upon a call by all others, receive z by broadcast from rank 0. */
{
   int size;

   MPI_Bcast (&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (size > 0) {
      _mpz_realloc (z, size);
      MPI_Bcast (z->_mp_d,  size, MPI_UNSIGNED_LONG, 0,
         MPI_COMM_WORLD);
   }
   else if (size < 0) {
      _mpz_realloc (z, -size);
      MPI_Bcast (z->_mp_d, -size, MPI_UNSIGNED_LONG, 0,
         MPI_COMM_WORLD);
   }
   z->_mp_size = size;
}

/*****************************************************************************/
/*                                                                           */
/* Submitting jobs and retrieving their results.                             */
/*                                                                           */
/*****************************************************************************/

void cm_mpi_broadcast_N (mpz_srcptr N)
   /* Send data depending on N to all workers. */
{
   int size, rank;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_BROADCAST_N,
         MPI_COMM_WORLD);
   mpi_bcast_send_mpz (N);
}

/*****************************************************************************/

void cm_mpi_broadcast_primorial (mpz_srcptr primorialB)
   /* Send the primorial to all workers. */
{
   int size, rank;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_BROADCAST_PRIMORIAL,
         MPI_COMM_WORLD);
   mpi_bcast_send_mpz (primorialB);
}

/*****************************************************************************/

void cm_mpi_broadcast_sqrt (int no_qstar, long int *qstar, mpz_t *qroot)
   /* Broadcast variables to all workers that are needed for the square
      root of discriminants and Cornacchia step.
      qstar and qroot are arrays of no_qstar signed primes and their
      square roots modulo N that the workers do not know yet; they append
      them to their list. */
{
   int size, rank;
   int i;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_BROADCAST_SQRT,
         MPI_COMM_WORLD);
   MPI_Bcast (&no_qstar, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (qstar, no_qstar, MPI_LONG, 0, MPI_COMM_WORLD);
   for (i = 0; i < no_qstar; i++)
      mpi_bcast_send_mpz (qroot [i]);
}

/*****************************************************************************/

void cm_mpi_clear_N ()
   /* Tell the workers to forget the incremental data sent for the square
      root step so that they can move on to the next N. */
{
   int size, rank;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_CLEAR_N,
         MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_submit_tonelli (int rank, int job, const long int a)
   /* Submit the Tonelli job of the given number to the worker of the given
      rank. The job number will be passed back by the worker so that the
      result can be identified. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_TONELLI, MPI_COMM_WORLD);
   MPI_Send (&a, 1, MPI_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_get_tonelli (mpz_ptr root, int rank, cm_stat_ptr stat)
   /* Get the result of a Tonelli job from worker rank and put it into root.
      Timing information from the worker is returned in stat. */
{
   mpi_recv_mpz (root, rank);
   MPI_Recv (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);
}

/*****************************************************************************/

void cm_mpi_submit_ecpp_one_step2 (int rank, int job, mpz_t *cert1,
   const char* modpoldir, bool tower)
   /* Submit the ECPP curve creation job of the given number to the worker
      of the given rank; the other parameters are as the input in
      cm_ecpp_one_step2 in ecpp.c */
{
   int tow, i;

   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_ECPP2, MPI_COMM_WORLD);
   for (i = 0; i < 4; i++)
      mpi_send_mpz (cert1 [i], rank);
   MPI_Send (modpoldir, strlen (modpoldir), MPI_CHAR, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD);
   tow = tower;
   MPI_Send (&tow, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_get_ecpp_one_step2 (mpz_t *cert2, int rank, cm_stat_ptr stat)
   /* Get the result of an ECPP curve creation job from worker rank and
      put it into cert2, as output by cm_ecpp_one_step2 in ecpp.c.
      Timing information from the worker is returned in stat. */
{
   int i;

   for (i = 0; i < 6; i++)
      mpi_recv_mpz (cert2 [i], rank);
   for (i = 1; i <= 3 ; i++)
   MPI_Recv (&(stat->timer [i]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);
}

/*****************************************************************************/

void cm_mpi_submit_curve_cardinalities (int rank, int job, int_cl_t *d,
   int no_d)
   /* Submit the job of the given number for determining a list of
      potential curve cardinalities for the array of no_d discriminants
      in d to the worker of the given rank. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_CARD, MPI_COMM_WORLD);
   MPI_Send (&no_d, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
   MPI_Send (d, no_d, MPI_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
}

/*****************************************************************************/

mpz_t* cm_mpi_get_curve_cardinalities (int *no_card, int_cl_t **card_d,
   int rank, cm_stat_ptr stat)
   /* Get the result of a curve cardinality job from worker rank.
      Timing information from the worker is returned in stat.
      Parameters and the return value are as in
      cm_ecpp_compute_cardinalities; in particular, card_d and the
      return value are newly allocated arrays. */
{
   mpz_t *res;
   int i;

   MPI_Recv (no_card, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
   *card_d = (int_cl_t *) malloc (*no_card * sizeof (int_cl_t));
   res = (mpz_t *) malloc (*no_card * sizeof (mpz_t));
   MPI_Recv (*card_d, *no_card, MPI_LONG, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);
   for (i = 0; i < *no_card; i++) {
      mpz_init (res [i]);
      mpi_recv_mpz (res [i], rank);
   }
   MPI_Recv (&(stat->timer [3]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);

   return res;
}

/*****************************************************************************/

void cm_mpi_submit_is_prime (int rank, int job, mpz_srcptr n)
   /* Submit the job of the given number for testing the primality of n
      to the worker of the given rank. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_PRIME, MPI_COMM_WORLD);
   mpi_send_mpz (n, rank);
}

/*****************************************************************************/

bool cm_mpi_get_is_prime (int rank, cm_stat_ptr stat)
   /* Get the result of a prime testing job from worker rank and return it.
      Timing information from the worker is returned in stat. */
{
   int res;

   MPI_Recv (&res, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
   MPI_Recv (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);

   return ((bool) res);
}

/*****************************************************************************/

void cm_mpi_submit_h_chunk (int rank, int job, uint_cl_t Dmin,
   uint_cl_t Dmax)
   /* Submit the job of the given number for computing a chunk of class
      numbers to the worker of the given rank; the other parameters are as
      the input of cm_ecpp_compute_h_chunk in ecpp.c. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_H, MPI_COMM_WORLD);
   MPI_Send (&Dmin, 1, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
   MPI_Send (&Dmax, 1, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_get_h_chunk (uint_cl_t *h, int rank, cm_stat_ptr stat)
   /* Get the result of a class number job from worker rank and
      put it into h, as output by cm_ecpp_compute_h_chunk in ecpp.c.
      Timing information from the worker is returned in stat. */
{
   MPI_Status status;
   int no;

   MPI_Probe (rank, MPI_TAG_DATA, MPI_COMM_WORLD, &status);
   MPI_Get_count (&status, MPI_UNSIGNED_LONG, &no);
   MPI_Recv (h, no, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD,
      NULL);
   MPI_Recv (&(stat->timer [5]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);
}

/*****************************************************************************/

void cm_mpi_submit_trial_div (int rank, int job, mpz_t *n, int no_n)
   /* Submit the job of the given number for trial dividing to the worker
      of the given rank; the other parameters are as the input of
      cm_ecpp_trial_div in ecpp.c, except that primorialB needs to have
      been sent independently. */
{
   int i;

   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_TRIAL_DIV, MPI_COMM_WORLD);
   MPI_Send (&no_n, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
   for (i = 0; i < no_n; i++)
      mpi_send_mpz (n [i], rank);
}

/*****************************************************************************/

void cm_mpi_get_trial_div (mpz_t *l, int rank, cm_stat_ptr stat)
   /* Get the result of a trial division job from worker rank and put it
      into l, as output by cm_ecpp_trial_div in ecpp.c.
      Timing information from the worker is returned in stat. */
{
   int no_n, i;

   MPI_Recv (&no_n, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
   for (i = 0; i < no_n; i++)
      mpi_recv_mpz (l [i], rank);
   MPI_Recv (&(stat->timer [1]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);
}

/*****************************************************************************/
/*                                                                           */
/* Worker implementation.                                                    */
/*                                                                           */
/*****************************************************************************/

static void mpi_worker ()
   /* The workers are started on rank 1 to size-1, and they run continually,
      accepting jobs until the server tells them to finish. */
{
   MPI_Status status;
   mpi_name_t name;
   int name_length;
   int job;
   bool finish;
   cm_stat_t stat;
   int i;
   
   /* Broadcast values. */
   mpz_t N, primorialB;
   int no_qstar, no_qstar_old, no_qstar_new;
   long int *qstar;
   mpz_t *qroot;

   /* Tonelli */
   long int a;
   unsigned int e = 0;
   mpz_t r, z, root;

   /* ECPP step 2 */
   mpz_t cert1 [4], cert2 [6];
   char *modpoldir;
   int len, tower;

   /* Curve cardinalities. */
   int_cl_t *d;
   int no_d, no_card;
   int_cl_t *card_d;
   mpz_t *card;

   /* Prime test. */
   mpz_t p;
   int isprime;

   /* Compute a chunk of h. */
   uint_cl_t Dmin, Dmax;
   uint_cl_t *h;
   int no;

   /* Trial division. */
   int no_n;
   mpz_t *l;

   mpz_init (N);
   mpz_init (primorialB);
   no_qstar = 0;
   qstar = (long int *) malloc (0);
   qroot = (mpz_t *) malloc (0);
   mpz_init (r);
   mpz_init (z);
   mpz_init (root);
   for (i = 0; i < 4; i++)
      mpz_init (cert1 [i]);
   for (i = 0; i < 6; i++)
      mpz_init (cert2 [i]);
   mpz_init (p);

   /* Gather data. */
   MPI_Get_processor_name (name, &name_length);
   MPI_Gather (name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
      NULL, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);

   finish = false;
   while (!finish) {
      /* Wait for a message from the server. */
      MPI_Recv (&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      cm_stat_init (stat);
      switch (status.MPI_TAG) {
      case MPI_TAG_JOB_BROADCAST_N:
         mpi_bcast_recv_mpz (N);
         e = cm_nt_mpz_tonelli_generator (r, z, N);
         break;
      case MPI_TAG_JOB_BROADCAST_PRIMORIAL:
         mpi_bcast_recv_mpz (primorialB);
         break;
      case MPI_TAG_JOB_BROADCAST_SQRT:
         MPI_Bcast (&no_qstar_new, 1, MPI_INT, 0, MPI_COMM_WORLD);
         no_qstar_old = no_qstar;
         no_qstar += no_qstar_new;
         qstar = (long int *)
            realloc (qstar, no_qstar * sizeof (long int));
         qroot = (mpz_t *) realloc (qroot, no_qstar * sizeof (mpz_t));
         MPI_Bcast (qstar + no_qstar_old, no_qstar_new, MPI_LONG, 0,
            MPI_COMM_WORLD);
         for (i = no_qstar_old; i < no_qstar; i++) {
            mpz_init (qroot [i]);
            mpi_bcast_recv_mpz (qroot [i]);
         }
         break;
      case MPI_TAG_JOB_CLEAR_N:
         for (i = 0; i < no_qstar; i++)
            mpz_clear (qroot [i]);
         no_qstar = 0;
         qstar = (long int *) realloc (qstar, 0);
         qroot = (mpz_t *) realloc (qroot, 0);
         break;
      case MPI_TAG_JOB_TONELLI:
         cm_timer_start (stat->timer [0]);
         /* Receive the input. */
         MPI_Recv (&a, 1, MPI_LONG, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);

         /* Compute the result. */
         cm_nt_mpz_tonelli_si_with_generator (root, a, N, e, r, z);

         /* Notify and send the result. */
         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_TONELLI, MPI_COMM_WORLD);
         mpi_send_mpz (root, 0);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         break;
      case MPI_TAG_JOB_ECPP2:
         for (i = 0; i < 4; i++)
            mpi_recv_mpz (cert1 [i], 0);
         MPI_Probe (0, MPI_TAG_DATA, MPI_COMM_WORLD, &status);
         MPI_Get_count (&status, MPI_CHAR, &len);
         modpoldir = (char *) malloc ((len + 1) * sizeof (char));
         MPI_Recv (modpoldir, len, MPI_CHAR, 0, MPI_TAG_DATA,
            MPI_COMM_WORLD, &status);
         modpoldir [len] = '\0';
         MPI_Recv (&tower, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);

         cm_ecpp_one_step2 (cert2, cert1, modpoldir, (bool) tower, true,
            false, stat);
         free (modpoldir);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_ECPP2, MPI_COMM_WORLD);
         for (i = 0; i < 6; i++)
            mpi_send_mpz (cert2 [i], 0);
         for (i = 1; i <= 3; i++)
            MPI_Send (&(stat->timer [i]->elapsed), 1, MPI_DOUBLE, 0,
               MPI_TAG_DATA, MPI_COMM_WORLD);
         break;
      case MPI_TAG_JOB_CARD:
         cm_timer_start (stat->timer [0]);
         MPI_Recv (&no_d, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            NULL);
         d = (int_cl_t *) malloc (no_d * sizeof (int_cl_t));
         MPI_Recv (d, no_d, MPI_LONG, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            NULL);

         card = cm_ecpp_compute_cardinalities (&no_card, &card_d, d, no_d,
            N, qstar, no_qstar, qroot);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_CARD, MPI_COMM_WORLD);
         MPI_Send (&no_card, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD);
         MPI_Send (card_d, no_card, MPI_LONG, 0, MPI_TAG_DATA,
            MPI_COMM_WORLD);
         for (i = 0; i < no_card; i++) {
            mpi_send_mpz (card [i], 0);
            mpz_clear (card [i]);
         }
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         free (d);
         free (card_d);
         free (card);
         break;
      case MPI_TAG_JOB_PRIME:
         cm_timer_start (stat->timer [0]);
         mpi_recv_mpz (p, 0);

         isprime = (int) cm_nt_is_prime (p);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_PRIME, MPI_COMM_WORLD);
         MPI_Send (&isprime, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         break;
      case MPI_TAG_JOB_H:
         cm_timer_start (stat->timer [0]);
         MPI_Recv (&Dmin, 1, MPI_UNSIGNED_LONG, 0, MPI_TAG_DATA,
            MPI_COMM_WORLD, &status);
         MPI_Recv (&Dmax, 1, MPI_UNSIGNED_LONG, 0, MPI_TAG_DATA,
            MPI_COMM_WORLD, &status);
         no = (Dmax - Dmin) / 2;
         h = (uint_cl_t *) malloc (no * sizeof (uint_cl_t));

         cm_ecpp_compute_h_chunk (h, Dmin, Dmax);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_H, MPI_COMM_WORLD);
         MPI_Send (h, no, MPI_UNSIGNED_LONG, 0, MPI_TAG_DATA,
            MPI_COMM_WORLD);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         free (h);
         break;
      case MPI_TAG_JOB_TRIAL_DIV:
         cm_timer_start (stat->timer [0]);
         MPI_Recv (&no_n, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);
         card = (mpz_t *) malloc (no_n * sizeof (mpz_t));
         l = (mpz_t *) malloc (no_n * sizeof (mpz_t));
         for (i = 0; i < no_n; i++) {
            mpz_init (card [i]);
            mpi_recv_mpz (card [i], 0);
            mpz_init (l [i]);
         }

         cm_ecpp_trial_div (l, card, no_n, primorialB);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_TRIAL_DIV,
            MPI_COMM_WORLD);
         MPI_Send (&no_n, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD);
         for (i = 0; i < no_n; i++)
            mpi_send_mpz (l [i], 0);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         for (i = 0; i < no_n; i++) {
            mpz_clear (card [i]);
            mpz_clear (l [i]);
         }
         free (card);
         free (l);
         break;
      case MPI_TAG_FINISH:
         finish = true;
         break;
      default:
         printf ("Unknown job type in mpi_worker; finishing\n");
         finish = true;
      }
   }

   mpz_clear (N);
   mpz_clear (primorialB);
   free (qstar);
   free (qroot);
   mpz_clear (r);
   mpz_clear (z);
   mpz_clear (root);
   mpz_clear (p);
   for (i = 0; i < 4; i++)
      mpz_clear (cert1 [i]);
   for (i = 0; i < 6; i++)
      mpz_clear (cert2 [i]);
}

/*****************************************************************************/
/*                                                                           */
/* Initialisations and worker queue handling.                                */
/*                                                                           */
/*****************************************************************************/

void cm_mpi_queue_push (int rank)
   /* Put worker of given rank back into the queue. */
{
   worker_queue [worker_queue_size] = rank;
   worker_queue_size++;
}

/*****************************************************************************/

int cm_mpi_queue_pop ()
   /* Get a worker rank from the queue, or return -1 if the queue
      is empty. */
{
   int rank;

   if (worker_queue_size == 0)
      rank = -1;
   else {
      worker_queue_size--;
      rank = worker_queue [worker_queue_size];
   }

   return rank;
}

/*****************************************************************************/

static void mpi_server_init (const int size, bool debug)
   /* The server is started on rank 0, initialises the workers and returns.
      The sequential code of the application should then be run in rank 0
      and occasionally make use of the workers for parallel sections. */
{
   mpi_name_t *worker_name;
   int name_length;
   int i;

   /* Set up worker queue. */
   worker_queue_size = size - 1;
   worker_queue = (int *) malloc (worker_queue_size * sizeof (int));
   for (i = 0; i < worker_queue_size; i++)
      worker_queue [i] = i+1;

   /* Gather worker names and set up queue of workers on local node. */
   worker_name = (mpi_name_t *) malloc (size * sizeof (mpi_name_t));
   MPI_Get_processor_name (worker_name [0], &name_length);
   MPI_Gather (MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
      worker_name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, 0, MPI_COMM_WORLD);
   worker_queue_local = (int *) malloc (worker_queue_size * sizeof (int));
   worker_queue_local_size = 0;
   for (i = 1; i < size; i++)
      if (!strncmp (worker_name [0], worker_name [i],
         MPI_MAX_PROCESSOR_NAME)) {
         worker_queue_local [worker_queue_local_size] = i;
         worker_queue_local_size++;
      }
   worker_queue_local = (int *) realloc (worker_queue_local,
      worker_queue_local_size * sizeof (int));
   free (worker_name);

   if (debug)
      printf ("MPI with %i workers initialised, of which %i are local.\n",
         worker_queue_size, worker_queue_local_size);
}

/*****************************************************************************/

static void mpi_server_clear (const int size)
{
   int i;
   const int dummy = 0;

   /* Tell the workers to finish. */
   for (i = 1; i < size; i++)
      MPI_Send (&dummy, 1, MPI_INT, i, MPI_TAG_FINISH, MPI_COMM_WORLD);

   free (worker_queue);
   free (worker_queue_local);
}

/*****************************************************************************/

void cm_mpi_init (bool debug)
   /* Start the MPI environment and launch the server and the workers. */
{
   int size, rank;

   MPI_Init (NULL, NULL);
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);

   if (rank == 0)
      mpi_server_init (size, debug);
   else
      mpi_worker ();
}

/*****************************************************************************/

void cm_mpi_clear ()
   /* Stop the server and the workers. */
{
   int size, rank;

   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if (rank == 0) {
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      mpi_server_clear (size);
   }

   MPI_Finalize ();
}

/*****************************************************************************/

