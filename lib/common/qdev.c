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
static long int lognormmax (ctype op);

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

static long int lognormmax (ctype op)
   /* computes the logarithm in base 2 (as the exponent of the value) of the */
   /* max norm of op                                                         */

{
   if (fsgn (op->re) == 0)
      if (fsgn (op->im) == 0)
         return fget_emin ();
      else
         return (fget_exp (op->im));
   else
      if (fsgn (op->im) == 0)
         return (fget_exp (op->re));
      else
         return MAX(fget_exp (op->re), fget_exp (op->im));
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_qdev_init (cm_qdev_t *f, fprec_t prec)
   /* initialises the addition chain for eta */

{
   int n, i, j, k;
//    int n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0, n6 = 0;

   f->length = 2 * ((fprec_t) (sqrt (prec * 0.085) + 1)) + 1;
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

void cm_qdev_eval (ctype rop, cm_qdev_t f, ctype q1)
   /* evaluates f in q */

{
   mp_exp_t prec;
   ctype    *q, term;
   int      n, i;

   prec = fget_prec (rop->re);

   q = (ctype *) malloc (f.length * sizeof (ctype));
   cinit (q [1], prec);
   cset (q [1], q1);
   cinit (term, prec);

   cset_si (rop, f.chain [0][5]);
   if (f.chain [1][5] != 0)
   {
      cmul_si (term, q [1], f.chain [1][5]);
      cadd (rop, rop, term);
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
         printf ("q "); cout_str (stdout, 10, 10, q [1]);
         printf ("\n");
         printf ("q^i "); cout_str (stdout, 10, 10, q [n-1]);
         printf ("\n");
         exit (1);
      }
      cinit (q [n], prec);
      switch (f.chain [n][1])
      {
      case 1:
         csqr (q [n], q [f.chain [n][2]]);
         break;
      case 2:
         cmul (q [n], q [f.chain [n][2]], q [f.chain [n][3]]);
         break;
      case 3:
         csqr (q [n], q [f.chain [n][2]]);
         csqr (q [n], q [n]);
         break;
      case 4:
         csqr (q [n], q [f.chain [n][2]]);
         cmul (q [n], q[n], q [f.chain [n][3]]);
         break;
      case 5:
         cmul (q [n], q [f.chain [n][2]], q [f.chain [n][3]]);
         csqr (q [n], q[n]);
         break;
      case 6:
         cmul (q [n], q [f.chain [n][2]], q [f.chain [n][3]]);
         cmul (q [n], q [n], q [f.chain [n][4]]);
         break;
      }
      if (f.chain [n][5] != 0)
      {
         cmul_si (term, q [n], f.chain [n][5]);
         cadd (rop, rop, term);
      }
   }

   for (i = 1; i <= n; i++)
      cclear (q [i]);
   free (q);
   cclear (term);
}

/*****************************************************************************/

void cm_qdev_eval_fr (ftype rop, cm_qdev_t f, ftype q1)
   /* evaluates f in q */

{
   mp_exp_t prec;
   ftype    *q, term;
   int      n, i;

   prec = fget_prec (rop);

   q = (ftype *) malloc (f.length * sizeof (ftype));
   finit (q [1], prec);
   fset (q [1], q1);
   finit (term, prec);

   fset_si (rop, f.chain [0][5]);
   fmul_si (term, q [1], f.chain [1][5]);
   fadd (rop, rop, term);
   n = 1;

   /* take next power of q into account if result is not precise enough */
   while (fget_exp (term) > -prec)
   {
      n++;
      if (n >= f.length)
      {
         printf ("*** Houston, we have a problem! Addition chain too short ");
         printf ("in 'qdev_eval_fr'.\n");
         printf ("n=%i, length=%i\n", n, f.length);
         printf ("q "); fout_str (stdout, 10, 0, q [1]);
         printf ("\n");
         printf ("q^i "); fout_str (stdout, 10, 0, q [n-1]);
         printf ("\n");
         exit (1);
      }
      finit (q [n], prec);
      switch (f.chain [n][1])
      {
      case 1:
         fsqr (q [n], q [f.chain [n][2]]);
         break;
      case 2:
         fmul (q [n], q [f.chain [n][2]], q [f.chain [n][3]]);
         break;
      case 3:
         fsqr (q [n], q [f.chain [n][2]]);
         fsqr (q [n], q [n]);
         break;
      case 4:
         fsqr (q [n], q [f.chain [n][2]]);
         fmul (q [n], q[n], q [f.chain [n][3]]);
         break;
      case 5:
         fmul (q [n], q [f.chain [n][2]], q [f.chain [n][3]]);
         fsqr (q [n], q[n]);
         break;
      case 6:
         fmul (q [n], q [f.chain [n][2]], q [f.chain [n][3]]);
         fmul (q [n], q [n], q [f.chain [n][4]]);
         break;
      }
      if (f.chain [n][5] != 0)
      {
         fmul_si (term, q [n], f.chain [n][5]);
         fadd (rop, rop, term);
      }
   }

   for (i = 1; i <= n; i++)
      fclear (q [i]);
   free (q);
   fclear (term);
}

/*****************************************************************************/
