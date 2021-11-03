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

#define MPI_TAG_READY       1
#define MPI_TAG_FINISH      2
#define MPI_TAG_DATA        3
#define MPI_TAG_JOB_TONELLI 4
#define MPI_TAG_JOB_ECPP2   5
#define MPI_TAG_JOB_CARD    6

static int *worker_queue;
static int worker_queue_size;

static void mpi_send_mpz (mpz_srcptr z, const int rank);
static void mpi_recv_mpz (mpz_ptr z, const int rank);
static void mpi_worker (const int rank);
static void mpi_server_init (const int size);
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
/*                                                                           */
/* Submitting jobs and retrieving their results.                             */
/*                                                                           */
/*****************************************************************************/

void cm_mpi_submit_tonelli (int rank, int job, const long int a,
   mpz_srcptr p, unsigned int e, mpz_srcptr q, mpz_srcptr z)
   /* Submit the Tonelli job of the given number to the worker of the given
      rank. The job number will be passed back by the worker so that the
      result can be identified. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_TONELLI, MPI_COMM_WORLD);
   MPI_Send (&a, 1, MPI_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
   mpi_send_mpz (p, rank);
   MPI_Send (&e, 1, MPI_UNSIGNED, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
   mpi_send_mpz (q, rank);
   mpi_send_mpz (z, rank);
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

void cm_mpi_submit_curve_cardinalities (int rank, int job, mpz_srcptr N,
   mpz_srcptr root, int_cl_t d)
   /* Submit the job of the given number for determining a list of 0 to 6
      potential curve cardinalities to the worker of the given rank; the
      other parameters are as the input in cm_ecpp_curve_cardinalities
      in ecpp.c. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_CARD, MPI_COMM_WORLD);
   mpi_send_mpz (N, rank);
   mpi_send_mpz (root, rank);
   MPI_Send (&d, 1, MPI_LONG, rank, MPI_TAG_DATA, MPI_COMM_WORLD);
}

/*****************************************************************************/

int cm_mpi_get_curve_cardinalities (mpz_t *n, int rank, cm_stat_ptr stat)
   /* Get the result of a curve cardinality job from worker rank and
      put it into n, as output by cm_ecpp_curve_cardinalities in ecpp.c.
      Timing information from the worker is returned in stat.
      The return value is the number of cardinalities. */
{
   int i, no;

   MPI_Recv (&no, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
   for (i = 0; i < no; i++)
      mpi_recv_mpz (n [i], rank);
   MPI_Recv (&(stat->timer [3]->elapsed), 1, MPI_DOUBLE, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD, NULL);

   return no;
}

/*****************************************************************************/
/*                                                                           */
/* Worker implementation.                                                    */
/*                                                                           */
/*****************************************************************************/

static void mpi_worker (const int rank)
   /* The workers are started on rank 1 to size-1, and they run continually,
      accepting jobs until the server tells them to finish. */
{
   MPI_Status status;
   int job;
   bool finish;
   cm_stat_t stat;
   int i;
   
   /* Tonelli */
   long int a;
   unsigned int e;
   mpz_t p, q, z, root;

   /* ECPP step 2 */
   mpz_t cert1 [4], cert2 [6];
   char *modpoldir;
   int len, tower;

   /* Curve cardinalities. */
   mpz_t N, n [6];
   int_cl_t d;
   int no;

   mpz_init (p);
   mpz_init (q);
   mpz_init (z);
   mpz_init (root);
   for (i = 0; i < 4; i++)
      mpz_init (cert1 [i]);
   for (i = 0; i < 6; i++) {
      mpz_init (cert2 [i]);
      mpz_init (n [i]);
   }
   mpz_init (N);
   mpz_init (root);

   /* Send a notification to the server. */
   MPI_Send (&rank, 1, MPI_INT, 0, MPI_TAG_READY, MPI_COMM_WORLD);

   finish = false;
   while (!finish) {
      /* Wait for a message from the server. */
      MPI_Recv (&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      cm_stat_init (stat);
      switch (status.MPI_TAG) {
      case MPI_TAG_JOB_TONELLI:
         cm_timer_start (stat->timer [0]);
         /* Receive the input. */
         MPI_Recv (&a, 1, MPI_LONG, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);
         mpi_recv_mpz (p, 0);
         MPI_Recv (&e, 1, MPI_UNSIGNED, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);
         mpi_recv_mpz (q, 0);
         mpi_recv_mpz (z, 0);

         /* Compute the result. */
         cm_nt_mpz_tonelli_si_with_generator (root, a, p, e, q, z);

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

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_TONELLI, MPI_COMM_WORLD);
         for (i = 0; i < 6; i++)
            mpi_send_mpz (cert2 [i], 0);
         for (i = 1; i <= 3; i++)
            MPI_Send (&(stat->timer [i]->elapsed), 1, MPI_DOUBLE, 0,
               MPI_TAG_DATA, MPI_COMM_WORLD);
         break;
      case MPI_TAG_JOB_CARD:
         cm_timer_start (stat->timer [0]);
         mpi_recv_mpz (N, 0);
         mpi_recv_mpz (root, 0);
         MPI_Recv (&d, 1, MPI_LONG, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);

         no = cm_ecpp_curve_cardinalities (n, N, root, d);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_CARD, MPI_COMM_WORLD);
         MPI_Send (&no, 1, MPI_INT, 0, MPI_TAG_DATA, MPI_COMM_WORLD);
         for (i = 0; i < no; i++)
            mpi_send_mpz (n [i], 0);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         break;
      case MPI_TAG_FINISH:
         finish = true;
         break;
      default:
         printf ("Unknown job type in mpi_worker; finishing\n");
         finish = true;
      }
   }

   mpz_clear (p);
   mpz_clear (q);
   mpz_clear (z);
   mpz_clear (root);
   for (i = 0; i < 4; i++)
      mpz_clear (cert1 [i]);
   for (i = 0; i < 6; i++) {
      mpz_clear (cert2 [i]);
      mpz_clear (n [i]);
   }
   mpz_clear (N);
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

static void mpi_server_init (const int size)
   /* The server is started on rank 0, initialises the workers and returns.
      The sequential code of the application should then be run in rank 0
      and occasionally make use of the workers for parallel sections. */
{
   int worker;

   /* Wait for all workers to join. */
   worker_queue = (int *) malloc ((size - 1) * sizeof (int));
   for (worker_queue_size = 0; worker_queue_size < size - 1;
      worker_queue_size++) {
      MPI_Recv (&worker, 1, MPI_INT, MPI_ANY_SOURCE, MPI_TAG_READY,
         MPI_COMM_WORLD, NULL);
      worker_queue [worker_queue_size] = worker;
   }
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
}

/*****************************************************************************/

void cm_mpi_init ()
   /* Start the MPI environment and launch the server and the workers. */
{
   int size, rank;

   MPI_Init (NULL, NULL);
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);

   if (rank == 0)
      mpi_server_init (size);
   else
      mpi_worker (rank);
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

