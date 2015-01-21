/*

qdev.c - code handling q-expansions

Copyright (C) 2009, 2015 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include "cm_common-impl.h"

#define MAX(a,b) ((a > b ? a : b))

static bool find_in_chain (int* index, cm_qdev_t f, int length, long int no);
static long int lognormmax (mpc_t op);

/*****************************************************************************/
/*                                                                           */
/* internal functions                                                        */
/*                                                                           */
/*****************************************************************************/

static bool find_in_chain (int* index, cm_qdev_t f, int length, long int no)
   /* looks for no in the first length elements of f                         */
   /* The return value indicates the success; if the operation succeeds,     */
   /* the index contains i such that f.chain [i][0] = no.                    */

{
   int left = 0, right = length - 1, middle;

   if (no < f.chain [0][0] || no > f.chain [length-1][0])
      return false;

   while (left < right - 1)
   {
      middle = (left + right) / 2;
      if (f.chain [middle][0] < no)
         left = middle;
      else
         right = middle;
   }
   if (f.chain [left][0] == no)
   {
      *index = left;
      return true;
   }
   else if (f.chain [right][0] == no)
   {
      *index = right;
      return true;
   }
   else
      return false;
}

/*****************************************************************************/

static long int lognormmax (mpc_t op)
   /* computes the logarithm in base 2 (as the exponent of the value) of the */
   /* max norm of op                                                         */

{
   if (mpfr_sgn (op->re) == 0)
      if (mpfr_sgn (op->im) == 0)
         return mpfr_get_emin ();
      else
         return (op->im->_mpfr_exp);
   else
      if (mpfr_sgn (op->im) == 0)
         return (op->re->_mpfr_exp);
      else
         return MAX(op->re->_mpfr_exp, op->im->_mpfr_exp);
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_qdev_init (cm_qdev_t *f, mp_prec_t prec)
   /* initialises the addition chain for eta */

{
   int n, i, j, k;
//    int n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0, n6 = 0;

   f->length = 2 * ((mp_prec_t) (sqrt (prec * 0.085) + 1)) + 1;
   /* must be odd                                                        */
   /* Since each power of q yields at least                              */
   /* log_2 (exp (sqrt (3) * pi)) = 7.85 bits,                           */
   /* the k yielding the largest exponent must satisfy                   */
   /* ((3*k+1)*k/2)*7.85 >= prec; we drop the +1 to simplify.            */
   /* Then we have twice as many exponents (taking into account the      */
   /* (3*k-1)*k/2), and one more for the constant coefficient.           */

   f->chain = (long int **) malloc (f->length * sizeof (long int *));
   for (n = 0; n < f->length; n++)
      f->chain [n] = (long int *) malloc (6 * sizeof (long int));

   (*f).chain [0][0] = 0;
   (*f).chain [0][5] = 1;
   for (n = 1; n <= f->length / 2; n++)
   {
      (*f).chain [2*n-1][0] = n*(3*n-1) / 2;
      (*f).chain [2*n][0] = n*(3*n+1) / 2;
      if (n % 2 == 0)
      {
         (*f).chain [2*n-1][5] = 1;
         (*f).chain [2*n][5] = 1;
      }
      else
      {
         (*f).chain [2*n-1][5] = -1;
         (*f).chain [2*n][5] = -1;
      }
   }

   (*f).chain [0][1] = 0;
   (*f).chain [1][1] = 0;
   for (n = 2; n < f->length; n++)
   {
      (*f).chain [n][1] = 0;
      /* try to express an even exponent as twice a previous one      */
      if ((*f).chain [n][0] % 2 == 0)
         if (find_in_chain (&i, *f, n, (*f).chain [n][0] / 2))
            {
               (*f).chain [n][1] = 1;
               (*f).chain [n][2] = i;
//                n1++;
            }
      /* try to express the exponent as the sum of two previous ones */
      for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
         if (find_in_chain (&j, *f, n, (*f).chain [n][0] - (*f).chain [i][0]))
            {
               (*f).chain [n][1] = 2;
               (*f).chain [n][2] = i;
               (*f).chain [n][3] = j;
//                n2++;
            }
      /* try to express an exponent as four times a previous one      */
      if ((*f).chain [n][0] % 4 == 0)
         if (find_in_chain (&i, *f, n, (*f).chain [n][0] / 4))
            {
               (*f).chain [n][1] = 3;
               (*f).chain [n][2] = i;
//                n3++;
            }
      /* try to express an even exponent as twice the sum of two previous */
      /* ones */
      /* try to express the exponent as twice a previous plus a third one */
      for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
         if (find_in_chain (&j, *f, n, (*f).chain [n][0] - 2 * (*f).chain [i][0]))
      {
         (*f).chain [n][1] = 4;
         (*f).chain [n][2] = i;
         (*f).chain [n][3] = j;
//          n4++;
      }
      if ((*f).chain [n][0] % 2 == 0)
         for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
            if (find_in_chain (&j, *f, n, (*f).chain [n][0]/2 - (*f).chain [i][0]))
               {
                  (*f).chain [n][1] = 5;
                  (*f).chain [n][2] = i;
                  (*f).chain [n][3] = j;
//                   n5++;
               }
      /* try to express the exponent as the sum of three previous ones */
      for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
         for (j = i; j < n && (*f).chain [n][1] == 0; j++)
            if (find_in_chain (&k, *f, n, (*f).chain [n][0] - (*f).chain [i][0]
                 - (*f).chain [j][0]))
            {
               (*f).chain [n][1] = 6;
               (*f).chain [n][2] = i;
               (*f).chain [n][3] = j;
               (*f).chain [n][4] = k;
//                n6++;
            }
      if ((*f).chain [n][1] == 0)
      {
         printf ("*** Houston, we have a problem! No success for element ");
         printf ("%i = %li in the addition chain ", n, (*f).chain [n][0]);
         printf ("computation in qdev_init. ");
         printf ("Go back programming!\n");
         exit (1);
      }
   }
/*
   printf ("n1 %i\n", n1);
   printf ("n2 %i\n", n2);
   printf ("n3 %i\n", n3);
   printf ("n4 %i\n", n4);
   printf ("n5 %i\n", n5);
   printf ("n6 %i\n", n6);
*/
}

/*****************************************************************************/

void cm_qdev_clear (cm_qdev_t *f)

{
   int n;

   for (n = 0; n < f->length; n++)
      free (f->chain [n]);
   free (f->chain);
}

/*****************************************************************************/

void cm_qdev_eval (mpc_t rop, cm_qdev_t f, mpc_t q1)
   /* evaluates f in q */

{
   mp_exp_t prec;
   mpc_t    *q, term;
   int      n, i;

   prec = mpfr_get_prec (rop->re);

   q = (mpc_t *) malloc (f.length * sizeof (mpc_t));
   mpc_init2 (q [1], prec);
   mpc_set (q [1], q1, MPC_RNDNN);
   mpc_init2 (term, prec);

   mpc_set_si (rop, f.chain [0][5], MPC_RNDNN);
   if (f.chain [1][5] != 0)
   {
      mpc_mul_si (term, q [1], f.chain [1][5], MPC_RNDNN);
      mpc_add (rop, rop, term, MPC_RNDNN);
   }
   n = 1;

   /* take next power of q into account if result is not precise enough */
   while (lognormmax (q [n]) > -prec)
   {
      n++;
      if (n >= f.length)
      {
         printf ("*** Houston, we have a problem! Addition chain too short ");
         printf ("in 'qdev_eval'.\n");
         printf ("n=%i, length=%i\n", n, f.length);
         printf ("q "); mpc_out_str (stdout, 10, 10, q [1], MPC_RNDNN);
         printf ("\n");
         printf ("q^i "); mpc_out_str (stdout, 10, 10, q [n-1], MPC_RNDNN);
         printf ("\n");
         exit (1);
      }
      mpc_init2 (q [n], prec);
      switch (f.chain [n][1])
      {
      case 1:
         mpc_sqr (q [n], q [f.chain [n][2]], MPC_RNDNN);
         break;
      case 2:
         mpc_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  MPC_RNDNN);
         break;
      case 3:
         mpc_sqr (q [n], q [f.chain [n][2]], MPC_RNDNN);
         mpc_sqr (q [n], q [n], MPC_RNDNN);
         break;
      case 4:
         mpc_sqr (q [n], q [f.chain [n][2]], MPC_RNDNN);
         mpc_mul (q [n], q[n], q [f.chain [n][3]], MPC_RNDNN);
         break;
      case 5:
         mpc_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  MPC_RNDNN);
         mpc_sqr (q [n], q[n], MPC_RNDNN);
         break;
      case 6:
         mpc_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  MPC_RNDNN);
         mpc_mul (q [n], q [n], q [f.chain [n][4]], MPC_RNDNN);
         break;
      }
      if (f.chain [n][5] != 0)
      {
         mpc_mul_si (term, q [n], f.chain [n][5], MPC_RNDNN);
         mpc_add (rop, rop, term, MPC_RNDNN);
      }
   }

   for (i = 1; i <= n; i++)
      mpc_clear (q [i]);
   free (q);
   mpc_clear (term);
}

/*****************************************************************************/

void cm_qdev_eval_fr (mpfr_t rop, cm_qdev_t f, mpfr_t q1)
   /* evaluates f in q */

{
   mp_exp_t prec;
   mpfr_t    *q, term;
   int      n, i;

   prec = mpfr_get_prec (rop);

   q = (mpfr_t *) malloc (f.length * sizeof (mpfr_t));
   mpfr_init2 (q [1], prec);
   mpfr_set (q [1], q1, GMP_RNDN);
   mpfr_init2 (term, prec);

   mpfr_set_si (rop, f.chain [0][5], GMP_RNDN);
   mpfr_mul_si (term, q [1], f.chain [1][5], GMP_RNDN);
   mpfr_add (rop, rop, term, GMP_RNDN);
   n = 1;

   /* take next power of q into account if result is not precise enough */
   while (mpfr_get_exp (term) > -prec)
   {
      n++;
      if (n >= f.length)
      {
         printf ("*** Houston, we have a problem! Addition chain too short ");
         printf ("in 'qdev_eval_fr'.\n");
         printf ("n=%i, length=%i\n", n, f.length);
         printf ("q "); mpfr_out_str (stdout, 10, 0, q [1], GMP_RNDN);
         printf ("\n");
         printf ("q^i "); mpfr_out_str (stdout, 10, 0, q [n-1], GMP_RNDN);
         printf ("\n");
         exit (1);
      }
      mpfr_init2 (q [n], prec);
      switch (f.chain [n][1])
      {
      case 1:
         mpfr_sqr (q [n], q [f.chain [n][2]], GMP_RNDN);
         break;
      case 2:
         mpfr_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  GMP_RNDN);
         break;
      case 3:
         mpfr_sqr (q [n], q [f.chain [n][2]], GMP_RNDN);
         mpfr_sqr (q [n], q [n], GMP_RNDN);
         break;
      case 4:
         mpfr_sqr (q [n], q [f.chain [n][2]], GMP_RNDN);
         mpfr_mul (q [n], q[n], q [f.chain [n][3]], GMP_RNDN);
         break;
      case 5:
         mpfr_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  GMP_RNDN);
         mpfr_sqr (q [n], q[n], GMP_RNDN);
         break;
      case 6:
         mpfr_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  GMP_RNDN);
         mpfr_mul (q [n], q [n], q [f.chain [n][4]], GMP_RNDN);
         break;
      }
      if (f.chain [n][5] != 0)
      {
         mpfr_mul_si (term, q [n], f.chain [n][5], GMP_RNDN);
         mpfr_add (rop, rop, term, GMP_RNDN);
      }
   }

   for (i = 1; i <= n; i++)
      mpfr_clear (q [i]);
   free (q);
   mpfr_clear (term);
}

/*****************************************************************************/
