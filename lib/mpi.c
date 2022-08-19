/*

mpicm.c - functions enabling MPI for CM

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

#define MPI_TAG_FINISH              1
#define MPI_TAG_DATA                2
#define MPI_TAG_JOB_PRIMORIAL       3
#define MPI_TAG_JOB_TONELLI         4
#define MPI_TAG_JOB_ECPP2           5
#define MPI_TAG_JOB_CARD            6
#define MPI_TAG_JOB_PRIME           7
#define MPI_TAG_JOB_H               8
#define MPI_TAG_JOB_TREE_GCD        9
#define MPI_TAG_JOB_BROADCAST_N    10
#define MPI_TAG_JOB_BROADCAST_SQRT 11
#define MPI_TAG_JOB_CLEAR_N        12
#define MPI_TAG_JOB_BROADCAST_INIT 13

typedef char mpi_name_t [MPI_MAX_PROCESSOR_NAME];

static int *worker_queue, *worker_queue_local;
static int worker_queue_size, worker_queue_local_size;
static MPI_Comm *split_comm;

static void mpi_send_mpz (mpz_srcptr z, const int rank,
   const MPI_Comm comm);
static void mpi_recv_mpz (mpz_ptr z, const int rank, const MPI_Comm comm);
static void mpi_bcast_send_mpz (mpz_srcptr z, const MPI_Comm comm);
static void mpi_bcast_recv_mpz (mpz_ptr z, const MPI_Comm comm);
static void mpi_worker (void);
static void mpi_server_init (const int size, bool debug);
static void mpi_server_clear (const int size);

static int mpi_compute_split_m (void);
static MPI_Comm mpi_communicator_split (void);
static void mpi_communicator_free (void);

/*****************************************************************************/
/*                                                                           */
/* Data structures sending and receiving.                                    */
/*                                                                           */
/*****************************************************************************/

static void mpi_send_mpz (mpz_srcptr z, const int rank, const MPI_Comm comm)
   /* Send z to rank in comm. */
{
   int size = z->_mp_size;

   MPI_Send (&size, 1, MPI_INT, rank, MPI_TAG_DATA, comm);
   if (size > 0)
      MPI_Send (z->_mp_d,  size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         comm);
   else if (size < 0)
      MPI_Send (z->_mp_d, -size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         comm);
}

/*****************************************************************************/

static void mpi_recv_mpz (mpz_ptr z, const int rank, const MPI_Comm comm)
   /* Get z from rank in comm. */
{
   int size;

   MPI_Recv (&size, 1, MPI_INT, rank, MPI_TAG_DATA, comm, NULL);
   if (size > 0) {
      _mpz_realloc (z, size);
      MPI_Recv (z->_mp_d,  size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         comm, NULL);
   }
   else if (size < 0) {
      _mpz_realloc (z, -size);
      MPI_Recv (z->_mp_d, -size, MPI_UNSIGNED_LONG, rank, MPI_TAG_DATA,
         comm, NULL);
   }
   z->_mp_size = size;
}

/*****************************************************************************/

static void mpi_bcast_send_mpz (mpz_srcptr z, MPI_Comm comm)
   /* Upon a call by rank 0, send z by broadcast to all others in comm. */
{
   int size = z->_mp_size;

   MPI_Bcast (&size, 1, MPI_INT, 0, comm);
   if (size > 0)
      MPI_Bcast (z->_mp_d,  size, MPI_UNSIGNED_LONG, 0, comm);
   else if (size < 0)
      MPI_Bcast (z->_mp_d, -size, MPI_UNSIGNED_LONG, 0, comm);
}

/*****************************************************************************/

static void mpi_bcast_recv_mpz (mpz_ptr z, MPI_Comm comm)
   /* Upon a call by all others in comm, receive z by broadcast from
      rank 0. */
{
   int size;

   MPI_Bcast (&size, 1, MPI_INT, 0, comm);
   if (size > 0) {
      _mpz_realloc (z, size);
      MPI_Bcast (z->_mp_d,  size, MPI_UNSIGNED_LONG, 0, comm);
   }
   else if (size < 0) {
      _mpz_realloc (z, -size);
      MPI_Bcast (z->_mp_d, -size, MPI_UNSIGNED_LONG, 0, comm);
   }
   z->_mp_size = size;
}

/*****************************************************************************/
/*                                                                           */
/* Helper functions for tree gcds.                                           */
/*                                                                           */
/*****************************************************************************/

static int mpi_compute_split_m ()
   /* Return the number m of communicators into which the world communicator
      is split; this is used for batching numbers in the trial division
      step. We would like to use about 16 communicators; since only
      m * floor (w / m) with w = size - 1  workers will be used for the gcd
      phase, it would also be nice if m divided w, but w + 1 is often a
      power of 2. So w try to choose m close to 16 so that w is not much
      above a multiple of m. */
{
   int size, w, i;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   w = size - 1;

   if (w < 16)
      return w;

   for (i = 16; i <= 18; i++)
      if (w % i == 0 || (w - 1) % i == 0 || (w - 2) % i == 0)
         return i;
   for (i = 15; i >= 14; i--)
      if (w % i == 0 || (w - 1) % i == 0 || (w - 2) % i == 0)
         return i;

   return 16;
}

/*****************************************************************************/

unsigned long int cm_mpi_compute_B ()
{
   int size, m;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   m = mpi_compute_split_m ();

   return (1ul << 29) * ((size - 1) / m);
}

/*****************************************************************************/

static MPI_Comm mpi_communicator_split ()
   /* Split the world communicator into m parts, as determined by
      mpi_compute_split_m, each of which contains (size - 1) / m
      consecutive workers, and put them into a global variable.
      Return the last communicator to which a process belongs;
      this is interesting only for the workers. */
{
   int m, size, w, group_size, i;
   MPI_Group world, group;
   int ranges [2][3];
   MPI_Comm res = MPI_COMM_NULL;

   m = mpi_compute_split_m ();
   split_comm = (MPI_Comm *) malloc (m * sizeof (MPI_Comm));

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   w = size - 1;
   group_size = w / m;

   MPI_Comm_group (MPI_COMM_WORLD, &world);
   ranges [0][0] = 0;
   ranges [0][1] = 0;
   ranges [0][2] = 1;
   ranges [1][2] = 1;
   for (i = 0; i < m; i++) {
      /* Put the processes 0 and i * group_size + 1, ..., (i+1) * group_size
         into communicator [i]. */
      ranges [1][0] = i * group_size + 1;
      ranges [1][1] = (i+1) * group_size;
      MPI_Group_range_incl (world, 2, ranges, &group);
      MPI_Comm_create (MPI_COMM_WORLD, group, split_comm + i);
      MPI_Group_free (&group);
      if (split_comm [i] != MPI_COMM_NULL)
         res = split_comm [i];
   }

   return res;
}

/*****************************************************************************/

static void mpi_communicator_free (void)
   /* Free the communicators created by mpi_communicator_split. */
{
   int m, i;

   m = mpi_compute_split_m ();

   for (i = 0; i < m; i++)
      if (split_comm [i] != MPI_COMM_NULL)
         MPI_Comm_free (split_comm + i);

   free (split_comm);
}


/*****************************************************************************/
/*                                                                           */
/* Submitting jobs and retrieving their results.                             */
/*                                                                           */
/*****************************************************************************/

void cm_mpi_broadcast_init (bool verbose, bool debug)
   /* Send common data to all workers. */
{
   int size, rank;
   int v, d;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_BROADCAST_INIT,
         MPI_COMM_WORLD);
   v = verbose;
   d = debug;
   MPI_Bcast (&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast (&d, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_broadcast_N (mpz_srcptr N)
   /* Send data depending on N to all workers. */
{
   int size, rank;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_BROADCAST_N,
         MPI_COMM_WORLD);
   mpi_bcast_send_mpz (N, MPI_COMM_WORLD);
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
      mpi_bcast_send_mpz (qroot [i], MPI_COMM_WORLD);
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

void cm_mpi_submit_primorial (char *tmpdir)
   /* Submit the job of computing primorials to all workers. */
{
   int size, rank;
   int len;

   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_PRIMORIAL,
         MPI_COMM_WORLD);
   len = (tmpdir ? strlen (tmpdir) + 1 : 0); /* +1 for trailing \0 */
   MPI_Bcast (&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (tmpdir)
      MPI_Bcast (tmpdir, len, MPI_CHAR, 0, MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_get_primorial (int rank, double *t)
   /* Get timing information of a primorial job from worker rank and return
      it in t. */
{
   MPI_Recv (t, 1, MPI_DOUBLE, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
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

void cm_mpi_get_tonelli (mpz_ptr root, int rank, double *t)
   /* Get the result of a Tonelli job from worker rank and put it into root.
      Timing information from the worker is returned in t. */
{
   mpi_recv_mpz (root, rank, MPI_COMM_WORLD);
   MPI_Recv (t, 1, MPI_DOUBLE, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
}

/*****************************************************************************/

void cm_mpi_submit_ecpp_one_step2 (int rank, int job, mpz_t *cert1,
   const char* modpoldir)
   /* Submit the ECPP curve creation job of the given number to the worker
      of the given rank; the other parameters are as the input in
      cm_ecpp_one_step2 in ecpp.c */
{
   int i;

   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_ECPP2, MPI_COMM_WORLD);
   for (i = 0; i < 4; i++)
      mpi_send_mpz (cert1 [i], rank, MPI_COMM_WORLD);
   MPI_Send (modpoldir, strlen (modpoldir), MPI_CHAR, rank, MPI_TAG_DATA,
      MPI_COMM_WORLD);
}

/*****************************************************************************/

void cm_mpi_get_ecpp_one_step2 (mpz_t *cert2, int rank, cm_stat_ptr stat)
   /* Get the result of an ECPP curve creation job from worker rank and
      put it into cert2, as output by cm_ecpp_one_step2 in ecpp.c.
      Timing information from the worker is returned in stat. */
{
   int i;

   for (i = 0; i < 6; i++)
      mpi_recv_mpz (cert2 [i], rank, MPI_COMM_WORLD);
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
   int rank, double *t)
   /* Get the result of a curve cardinality job from worker rank.
      Timing information from the worker is returned in t.
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
      mpi_recv_mpz (res [i], rank, MPI_COMM_WORLD);
   }
   MPI_Recv (t, 1, MPI_DOUBLE, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);

   return res;
}

/*****************************************************************************/

void cm_mpi_submit_is_prime (int rank, int job, mpz_srcptr n)
   /* Submit the job of the given number for testing the primality of n
      to the worker of the given rank. */
{
   MPI_Send (&job, 1, MPI_INT, rank, MPI_TAG_JOB_PRIME, MPI_COMM_WORLD);
   mpi_send_mpz (n, rank, MPI_COMM_WORLD);
}

/*****************************************************************************/

bool cm_mpi_get_is_prime (int rank, double *t)
   /* Get the result of a prime testing job from worker rank and return it.
      Timing information from the worker is returned in t. */
{
   int res;

   MPI_Recv (&res, 1, MPI_INT, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
   MPI_Recv (t, 1, MPI_DOUBLE, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);

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

void cm_mpi_get_h_chunk (unsigned int *h, int rank, double *t)
   /* Get the result of a class number job from worker rank and
      put it into h, as output by cm_ecpp_compute_h_chunk in ecpp.c.
      Timing information from the worker is returned in t. */
{
   MPI_Status status;
   int no;

   MPI_Probe (rank, MPI_TAG_DATA, MPI_COMM_WORLD, &status);
   MPI_Get_count (&status, MPI_UNSIGNED, &no);
   MPI_Recv (h, no, MPI_UNSIGNED, rank, MPI_TAG_DATA, MPI_COMM_WORLD,
      NULL);
   MPI_Recv (t, 1, MPI_DOUBLE, rank, MPI_TAG_DATA, MPI_COMM_WORLD, NULL);
}

/*****************************************************************************/

void cm_mpi_submit_tree_gcd (mpz_t *n, int no_n)
   /* Submit to all workers the jobs of executing cm_mpz_tree_gcd with the
      known primorial and the remaining arguments. */
{
   int size, rank;
   MPI_Comm comm;
   int m, offset, len, i, j;

   /* Alert all workers of the phase, since they are waiting
      for a message in MPI_COMM_WORLD. */
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   for (rank = 1; rank < size; rank++)
      MPI_Send (&rank, 1, MPI_INT, rank, MPI_TAG_JOB_TREE_GCD,
      MPI_COMM_WORLD);

   /* Split the numbers to be trial divided into m batches and
      send them to the different communicators. */
   m = mpi_compute_split_m ();
   len = no_n / m;
   for (i = 0; i < m; i++) {
      comm = split_comm [i];
      offset = i * len;
      if (i == m - 1)
         len = no_n - offset;
      MPI_Bcast (&len, 1, MPI_INT, 0, comm);
      for (j = 0; j < len; j++)
         mpi_bcast_send_mpz (n [offset + j], comm);
   }
}

/*****************************************************************************/

void cm_mpi_get_tree_gcd (mpz_t *gcd, int no_n, double *t)
   /* Get from all workers the results of gcd jobs as output by
      cm_mpz_tree_gcd and collect them into gcd. */
{
   int size, rank, job;
   MPI_Comm comm;
   MPI_Status status;
   mpz_t tmp;
   double t_local;
   int m, offset, len, i, j, k;

   *t = 0;
   mpz_init (tmp);

   m = mpi_compute_split_m ();
   for (i = 0; i < m; i++) {
      /* Handle batch i of numbers. */
      comm = split_comm [i];
      MPI_Comm_size (comm, &size);
      len = no_n / m;
      offset = i * len;
      if (i == m - 1)
         len = no_n - offset;
      /* Different workers compute gcds with different divisors of the
         primorial; we need to multiply them together. */
      for (k = 0; k < len; k++)
         mpz_set_ui (gcd [offset + k], 1);
      for (j = 1; j < size; j++) {
         MPI_Recv (&job, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
            comm, &status);
         rank = status.MPI_SOURCE;
         for (k = 0; k < len; k++) {
            mpi_recv_mpz (tmp, rank, comm);
            mpz_mul (gcd [offset + k], gcd [offset + k], tmp);
         }
         MPI_Recv (&t_local, 1, MPI_DOUBLE, rank, MPI_TAG_DATA, comm, NULL);
         *t += t_local;
      }
   }
   mpz_clear (tmp);
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
   int size, rank, tree_rank, job, flag;
   MPI_Comm tree_comm;
   bool finish;
   cm_stat_t stat;
   int i;
   
   /* Broadcast values. */
   int b;
   bool verbose = false, debug = false;
   mpz_t N;
   int no_qstar, no_qstar_old, no_qstar_new;
   long int *qstar;
   mpz_t *qroot;

   /* Primorial. */
   int len;
   char *tmpdir;
   mpz_t prim;
   bool read;

   /* Tonelli */
   long int a;
   unsigned int e = 0;
   mpz_t r, z, root;

   /* ECPP step 2 */
   mpz_t cert1 [4], cert2 [6];
   char *modpoldir;

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
   unsigned int *h;
   int no;

   /* Tree gcd. */
   int no_m;
   mpz_t *m, *gcd;

   mpz_init (N);
   mpz_init (prim);
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

   /* Split communicator. */
   tree_comm = mpi_communicator_split ();
   MPI_Comm_size (MPI_COMM_WORLD, &size);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if (tree_comm != MPI_COMM_NULL)
      MPI_Comm_rank (tree_comm, &tree_rank);
   else
      tree_rank = -1;

   /* Gather data. */
   MPI_Get_processor_name (name, &name_length);
   MPI_Gather (name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
      NULL, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);

   finish = false;
   while (!finish) {
      /* Wait for a message from the server, but avoid being busy by
         adding microsleep. */
      do {
         usleep (100);
         MPI_Iprobe (0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag,
         MPI_STATUS_IGNORE);
      } while (!flag);
      MPI_Recv (&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      cm_stat_init (stat);
      switch (status.MPI_TAG) {
      case MPI_TAG_JOB_BROADCAST_INIT:
         MPI_Bcast (&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
         verbose = b;
         MPI_Bcast (&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
         debug = b;
         break;
      case MPI_TAG_JOB_BROADCAST_N:
         mpi_bcast_recv_mpz (N, MPI_COMM_WORLD);
         e = cm_nt_mpz_tonelli_generator (r, z, N);
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
            mpi_bcast_recv_mpz (qroot [i], MPI_COMM_WORLD);
         }
         break;
      case MPI_TAG_JOB_CLEAR_N:
         for (i = 0; i < no_qstar; i++)
            mpz_clear (qroot [i]);
         no_qstar = 0;
         qstar = (long int *) realloc (qstar, 0);
         qroot = (mpz_t *) realloc (qroot, 0);
         break;
      case MPI_TAG_JOB_PRIMORIAL:
         cm_timer_start (stat->timer [0]);
         MPI_Bcast (&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
         if (len > 0) {
            tmpdir = (char *) malloc (len * sizeof (char));
            MPI_Bcast (tmpdir, len, MPI_CHAR, 0, MPI_COMM_WORLD);
         }
         else
            tmpdir = 0;

         /* The numbers to be trial divided are split into m batches,
            each of which is handled by the workers behind one of the
            split_comm. If split_size denotes the size of the communicator
            exluding rank 0 (that is, split_size = floor ((size - 1) / m)),
            then each worker inside chooses an interval of 2^29 and thus
            handles a product of value about exp (2^29), or 100MB.
            Currently this implies that the effective value of B is
            split_size * 2^29 = floor ((size - 1) / m) * 2^29.
            If size - 1 is not divisible by m, the last few workers
            do not take part in the computation. */
         if (tree_rank >= 1) {
            if (tmpdir)
               read = cm_file_read_primorial (tmpdir, prim, tree_rank - 1);
            else
               read = false;
            if (!read)
               cm_pari_prime_product (prim,
                  ((unsigned long int) tree_rank - 1) << 29,
                  ((unsigned long int) tree_rank)     << 29);
            if (tmpdir && !read && rank == tree_rank)
               /* Let only the workers in the first communicator write
                  the primorial to disk. */
               cm_file_write_primorial (tmpdir, prim, rank - 1);
         }

         if (len > 0)
            free (tmpdir);
         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_PRIMORIAL,
            MPI_COMM_WORLD);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
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
         mpi_send_mpz (root, 0, MPI_COMM_WORLD);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         break;
      case MPI_TAG_JOB_ECPP2:
         for (i = 0; i < 4; i++)
            mpi_recv_mpz (cert1 [i], 0, MPI_COMM_WORLD);
         MPI_Probe (0, MPI_TAG_DATA, MPI_COMM_WORLD, &status);
         MPI_Get_count (&status, MPI_CHAR, &len);
         modpoldir = (char *) malloc ((len + 1) * sizeof (char));
         MPI_Recv (modpoldir, len, MPI_CHAR, 0, MPI_TAG_DATA,
            MPI_COMM_WORLD, &status);
         modpoldir [len] = '\0';

         cm_ecpp_one_step2 (cert2, cert1, job, modpoldir, verbose, debug,
            stat);
         free (modpoldir);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_ECPP2, MPI_COMM_WORLD);
         for (i = 0; i < 6; i++)
            mpi_send_mpz (cert2 [i], 0, MPI_COMM_WORLD);
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
            mpi_send_mpz (card [i], 0, MPI_COMM_WORLD);
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
         mpi_recv_mpz (p, 0, MPI_COMM_WORLD);

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
         h = (unsigned int *) malloc (no * sizeof (unsigned int));

         cm_ecpp_compute_h_chunk (h, Dmin, Dmax);

         MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_H, MPI_COMM_WORLD);
         MPI_Send (h, no, MPI_UNSIGNED, 0, MPI_TAG_DATA, MPI_COMM_WORLD);
         cm_timer_stop (stat->timer [0]);
         MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
            MPI_TAG_DATA, MPI_COMM_WORLD);
         free (h);
         break;
      case MPI_TAG_JOB_TREE_GCD:
         if (tree_comm != MPI_COMM_NULL) {
            cm_timer_start (stat->timer [0]);
            MPI_Bcast (&no_m, 1, MPI_INT, 0, tree_comm);
            m = (mpz_t *) malloc (no_m * sizeof (mpz_t));
            gcd = (mpz_t *) malloc (no_m * sizeof (mpz_t));
            for (i = 0; i < no_m; i++) {
               mpz_init (m [i]);
               mpi_bcast_recv_mpz (m [i], tree_comm);
               mpz_init (gcd [i]);
            }

            cm_nt_mpz_tree_gcd (gcd, prim, m, no_m);

            MPI_Send (&job, 1, MPI_INT, 0, MPI_TAG_JOB_TREE_GCD, tree_comm);
            for (i = 0; i < no_m; i++) {
               mpi_send_mpz (gcd [i], 0, tree_comm);
               mpz_clear (gcd [i]);
               mpz_clear (m [i]);
            }

            cm_timer_stop (stat->timer [0]);
            MPI_Send (&(stat->timer [0]->elapsed), 1, MPI_DOUBLE, 0,
               MPI_TAG_DATA, tree_comm);
            free (m);
            free (gcd);
         }
         break;
      case MPI_TAG_FINISH:
         finish = true;
         break;
      default:
         printf ("Unknown job type in mpi_worker; finishing\n");
         finish = true;
      }
   }

   /* Free split communicators. */
   mpi_communicator_free ();

   mpz_clear (N);
   mpz_clear (prim);
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
   /* The server is started on rank 0, initialises the workers and
      communicators and returns.
      The sequential code of the application should then be run in rank 0
      and occasionally make use of the workers for parallel sections. */
{
   mpi_name_t *worker_name;
   int name_length;
   int i;

   /* Set up the split communicators. */
   mpi_communicator_split ();

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
      printf ("MPI with %i workers in %i groups initialised; "
              "%i workers are local.\n",
         worker_queue_size, mpi_compute_split_m (),
         worker_queue_local_size);
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

   mpi_communicator_free ();
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

