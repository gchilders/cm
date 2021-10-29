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

#define mpi_log(rk, ...) {printf ("[%2i] ", rk); printf (__VA_ARGS__);}

static int *worker_queue;
static int worker_queue_size;

static void mpi_send_mpz (mpz_srcptr z, const int rank);
static void mpi_recv_mpz (mpz_ptr z, const int rank);
static void mpi_worker (const int rank, bool debug);
static void mpi_server_init (const int size, bool debug);
static void mpi_server_clear (const int size, bool debug);

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

void cm_mpi_get_tonelli (mpz_ptr root, int rank)
   /* Get the result of a Tonelli job from worker rank and put it
      into root. */
{
   mpi_recv_mpz (root, rank);
}

/*****************************************************************************/
/*                                                                           */
/* Worker implementation.                                                    */
/*                                                                           */
/*****************************************************************************/

static void mpi_worker (const int rank, bool debug)
   /* The workers are started on rank 1 to size-1, and they run continually,
      accepting jobs until the server tells them to finish. */
{
   MPI_Status status;
   int job;
   bool finish;
   
   /* Tonelli */
   long int a;
   unsigned int e;
   mpz_t p, q, z, root;

   mpz_init (p);
   mpz_init (q);
   mpz_init (z);
   mpz_init (root);

   /* Send a notification to the server. */
   if (debug)
      mpi_log (rank, "started\n");
   MPI_Send (&rank, 1, MPI_INT, 0, MPI_TAG_READY, MPI_COMM_WORLD);

   finish = false;
   while (!finish) {
      /* Wait for a message from the server. */
      MPI_Recv (&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      switch (status.MPI_TAG) {
      case MPI_TAG_JOB_TONELLI:
         /* Receive the input. */
         MPI_Recv (&a, 1, MPI_LONG, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);
         mpi_recv_mpz (p, 0);
         MPI_Recv (&e, 1, MPI_UNSIGNED, 0, MPI_TAG_DATA, MPI_COMM_WORLD,
            &status);
         mpi_recv_mpz (q, 0);
         mpi_recv_mpz (z, 0);

         /* Compute the result. */
         if (debug)
            mpi_log (rank, "Tonelli %i\n", job);
         cm_nt_mpz_tonelli_si_with_generator (root, a, p, e, q, z);

         /* Notify and send the result. */
         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_TONELLI, MPI_COMM_WORLD);
         mpi_send_mpz (root, 0);
         break;
      case MPI_TAG_FINISH:
         if (debug)
            mpi_log (rank, "finished\n");
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
   int worker;

   /* Wait for all workers to join. */
   if (debug)
      mpi_log (0, "started\n");
   worker_queue = (int *) malloc ((size - 1) * sizeof (int));
   for (worker_queue_size = 0; worker_queue_size < size - 1;
      worker_queue_size++) {
      MPI_Recv (&worker, 1, MPI_INT, MPI_ANY_SOURCE, MPI_TAG_READY,
         MPI_COMM_WORLD, NULL);
      worker_queue [worker_queue_size] = worker;
      if (debug)
         mpi_log (0, "%i added to queue\n", worker);
   }
}

/*****************************************************************************/

static void mpi_server_clear (const int size, bool debug)
{
   int i;
   const int dummy = 0;

   /* Tell the workers to finish. */
   for (i = 1; i < size; i++) {
      MPI_Send (&dummy, 1, MPI_INT, i, MPI_TAG_FINISH, MPI_COMM_WORLD);
      if (debug)
         mpi_log (0, "%i removed from queue\n", i);
   }
   if (debug)
      mpi_log (0, "finished\n");

   free (worker_queue);
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
      mpi_worker (rank, debug);
}

/*****************************************************************************/

void cm_mpi_clear (bool debug)
   /* Stop the server and the workers. */
{
   int size, rank;

   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if (rank == 0) {
      MPI_Comm_size (MPI_COMM_WORLD, &size);
      mpi_server_clear (size, debug);
   }

   MPI_Finalize ();
}

/*****************************************************************************/

