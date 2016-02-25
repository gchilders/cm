/*

qdev.c - code handling q-expansions

Copyright (C) 2009, 2015, 2016 Andreas Enge

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

static bool find_in_chain (int* index, cm_qdev_t f, int length, long int no);
static double lognorm2 (ctype op);

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

static double lognorm2 (ctype op)
   /* computes the logarithm in base 2 of the complex norm of op */

{
   double   re, im;
   long int ere, eim, diff;

   /* Just extracting a double may overflow, so treat the exponents
      separately. */
   re = fget_d_2exp (&ere, crealref (op));
   im = fget_d_2exp (&eim, cimagref (op));

   /* Handle the case of 0 in one part separately, as it may be coupled
      with another very small exponent; then normalising for the larger
      exponent yields 0 and a problem with the logarithm. */
   if (re == 0)
      return (eim + log2 (fabs (im)));
   else if (im == 0)
      return (ere + log2 (fabs (re)));

   /* Normalise to keep the larger exponent; the smaller one may underflow,
      then the number becomes a harmless 0. */
   if (ere > eim) {
      diff = ere - eim;
      eim = ere;
      im /= (1ul << diff);
   }
   else {
      diff = eim - ere;
      ere = eim;
      re /= (1ul << diff);
   }

   return (ere + log2 (re*re + im*im) / 2);
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_qdev_init (cm_qdev_t *f, fprec_t prec)
   /* initialises the addition chain for eta */

{
   int n, i, j;

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
      f->chain [n] = (long int *) malloc (5 * sizeof (long int));

   f->chain [0][0] = 0;
   f->chain [0][4] = 1;
   for (n = 1; n <= f->length / 2; n++)
   {
      f->chain [2*n-1][0] = n*(3*n-1) / 2;
      f->chain [2*n][0] = n*(3*n+1) / 2;
      if (n % 2 == 0)
      {
         f->chain [2*n-1][4] = 1;
         f->chain [2*n][4] = 1;
      }
      else
      {
         f->chain [2*n-1][4] = -1;
         f->chain [2*n][4] = -1;
      }
   }

   f->chain [0][1] = 0;
   f->chain [1][1] = 0;
   for (n = 2; n < f->length; n++)
   {
      f->chain [n][1] = 0;
      /* try to express an even exponent as twice a previous one      */
      if (f->chain [n][0] % 2 == 0)
         if (find_in_chain (&i, *f, n, f->chain [n][0] / 2))
            {
               f->chain [n][1] = 1;
               f->chain [n][2] = i;
            }
      /* try to express the exponent as the sum of two previous ones */
      for (i = 0; i < n && f->chain [n][1] == 0; i++)
         if (find_in_chain (&j, *f, n, f->chain [n][0] - f->chain [i][0]))
            {
               f->chain [n][1] = 2;
               f->chain [n][2] = i;
               f->chain [n][3] = j;
            }
      /* try to express the exponent as twice a previous plus a third one */
      for (i = 0; i < n && f->chain [n][1] == 0; i++)
         if (find_in_chain (&j, *f, n, f->chain [n][0] - 2 * f->chain [i][0]))
      {
         f->chain [n][1] = 3;
         f->chain [n][2] = i;
         f->chain [n][3] = j;
      }
      /* This covers all cases for eta, see Enge-Johansson 2016. */
   }
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
   /* evaluates f in q1 */

{
   mp_prec_t prec;
   long int  local_prec;
   double    delta;
   ctype     *q, term, tmp1, tmp2;
   int       n, i;

   prec = fget_prec (rop->re);
   delta = - lognorm2 (q1);

   q = (ctype *) malloc (f.length * sizeof (ctype));
   cinit (q [1], prec);
   cset (q [1], q1);
   cinit (term, prec);
   cinit (tmp1, prec);
   cinit (tmp2, prec);

   cset_si (rop, f.chain [0][4]);
   if (f.chain [1][4] == 1)
     cadd (rop, rop, q [1]);
   else if (f.chain [1][4] == -1)
     csub (rop, rop, q [1]);
   else if (f.chain [1][4] != 0)
   {
      cmul_si (term, q [1], f.chain [1][4]);
      cadd (rop, rop, term);
   }

   n = 2;
   /* Adapt the precision for the next term. */
   local_prec = (long int) prec - (long int) (f.chain [n][0] * delta);

   while (local_prec >= 2)
   {
      cinit (q [n], (mp_prec_t) local_prec);
      switch (f.chain [n][1])
      {
      case 1:
         /* Reduce the precision of the argument to save some more time. */
         cset_prec (tmp1, local_prec);
         cset (tmp1, q [f.chain [n][2]]);
         csqr (q [n], tmp1);
         break;
      case 2:
         cset_prec (tmp1, local_prec);
         cset_prec (tmp2, local_prec);
         cset (tmp1, q [f.chain [n][2]]);
         cset (tmp2, q [f.chain [n][3]]);
         cmul (q [n], tmp1, tmp2);
         break;
      case 3:
         cset_prec (tmp1, local_prec);
         cset_prec (tmp2, local_prec);
         cset (tmp1, q [f.chain [n][2]]);
         cset (tmp2, q [f.chain [n][3]]);
         csqr (q [n], tmp1);
         cmul (q [n], q [n], tmp2);
         break;
      }
      if (f.chain [n][4] == 1)
        cadd (rop, rop, q [n]);
      else if (f.chain [n][4] == -1)
        csub (rop, rop, q [n]);
      else if (f.chain [n][4] != 0)
      {
	 cset_prec (term, (mp_prec_t) local_prec);
         cmul_si (term, q [n], f.chain [n][4]);
         cadd (rop, rop, term);
      }
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
      local_prec = (long int) prec - (long int) (f.chain [n][0] * delta);
   }

   for (i = 1; i < n; i++)
      cclear (q [i]);
   free (q);
   cclear (term);
   cclear (tmp1);
   cclear (tmp2);
}

/*****************************************************************************/

void cm_qdev_eval_fr (ftype rop, cm_qdev_t f, ftype q1)
   /* evaluates f in q1 */

{
   mp_prec_t prec;
   long int  local_prec, e;
   double    mantissa, delta;
   ftype     *q, term;
   int       n, i;

   prec = fget_prec (rop);
   mantissa = fget_d_2exp (&e, q1);
   delta = - (e + log2 (fabs (mantissa)));

   q = (ftype *) malloc (f.length * sizeof (ftype));
   finit (q [1], prec);
   fset (q [1], q1);
   finit (term, prec);

   fset_si (rop, f.chain [0][4]);
   if (f.chain [1][4] == 1)
     fadd (rop, rop, q [1]);
   else if (f.chain [1][4] == -1)
     fsub (rop, rop, q [1]);
   else if (f.chain [1][4] != 0)
   {
      fmul_si (term, q [1], f.chain [1][4]);
      fadd (rop, rop, term);
   }

   n = 2;
   /* Adapt the precision for the next term. */
   local_prec = (long int) prec - (long int) (f.chain [n][0] * delta);

   while (local_prec >= 2)
   {
      finit (q [n], (mp_prec_t) local_prec);
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
         fmul (q [n], q[n], q [f.chain [n][3]]);
         break;
      }
      if (f.chain [n][4] == 1)
        fadd (rop, rop, q [n]);
      else if (f.chain [n][4] == -1)
        fsub (rop, rop, q [n]);
      else if (f.chain [n][4] != 0)
      {
	 fset_prec (term, (mp_prec_t) local_prec);
         fmul_si (term, q [n], f.chain [n][4]);
         fadd (rop, rop, term);
      }
      n++;
      if (n >= f.length)
      {
         printf ("*** Houston, we have a problem! Addition chain too short ");
         printf ("in 'qdev_eval_fr'.\n");
         printf ("n=%i, length=%i\n", n, f.length);
         printf ("q "); fout_str (stdout, 10, 10, q [1]);
         printf ("\n");
         printf ("q^i "); fout_str (stdout, 10, 10, q [n-1]);
         printf ("\n");
         exit (1);
      }
      local_prec = (long int) prec - (long int) (f.chain [n][0] * delta);
   }

   for (i = 1; i < n; i++)
      fclear (q [i]);
   free (q);
   fclear (term);
}

/*****************************************************************************/

