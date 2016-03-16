/*

qdev.c - code handling q-expansions

Copyright (C) 2009, 2015, 2016 Andreas Enge

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

#include "cm_common-impl.h"

static bool find_in_chain (int* index, cm_qdev_t f, int length, long int no);
static double lognorm2 (ctype op);
static void qdev_eval_addition_sequence (ctype rop, cm_qdev_t f, ctype q1,
   double delta, int N);
static int dense_addition_sequence (long int ** b, int m);
static int minimal_dense_addition_sequence (long int ** b, int m);
static void qdev_eval_bsgs (ctype rop, cm_qdev_t f, ctype q1,
   double delta, int N);

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

static void qdev_eval_addition_sequence (ctype rop, cm_qdev_t f, ctype q1,
   double delta, int N)
   /* Evaluate f in q1 using the optimised addition sequence from f.
      N is the last index used in the addition chain.
      delta is the number of bits gained with each power of q.
      rop and q1 may be the same. */

{
   mp_prec_t prec;
   long int  local_prec;
   ctype     *q, term, tmp1, tmp2;
   int       n, i;

   prec = fget_prec (crealref (rop));

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

   for (n = 2; n <= N; n++) {
      local_prec = (long int) prec - (long int) (f.chain [n][0] * delta);
      cinit (q [n], (mp_prec_t) local_prec);
      switch (f.chain [n][1]) {
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
      else if (f.chain [n][4] != 0) {
	 cset_prec (term, (mp_prec_t) local_prec);
         cmul_si (term, q [n], f.chain [n][4]);
         cadd (rop, rop, term);
      }
   }

   for (i = 1; i < n; i++)
      cclear (q [i]);
   free (q);
   cclear (term);
   cclear (tmp1);
   cclear (tmp2);
}

/*****************************************************************************/

static int dense_addition_sequence (long int ** b, int m)
   /* Compute an addition sequence for b and return it via b itself.
      The return value is the number of additional entries.
      b is a matrix of m+1 rows and (at least) 2 columns; the row i
      represents the exponent i. Initially, b [i][0] is expected to be 1
      if the exponent i occurs, 0 otherwise
      During the course of the algorithm, new elements are added;
      these are marked by a 2 in b [i][0].
      b [i][1] is set to j if the exponent occurs and can be written as
        j+(i-j) with also occurring exponents j and i-j.
      The algorithm prefers doublings over other additions. */
{
   int i, j, found, added;

   added = 0;
   for (i = m; i >= 2; i--)
      if (b [i][0] != 0) {
         found = 0;
         /* Search for two existing entries adding to i. */
         for (j = i / 2; j >= 1 && !found; j--)
            /* In this way, a doubling is found if it exists. */
            if (b [j][0] != 0 && b [i-j][0] != 0) {
               b [i][1] = j;
               found = 1;
            }
         if (!found) {
            /* Add missing elements. The middle strategy seems to give the
               best performance; when minimising additionally, all three
               have a very similar outcome. */
#if 0
            /* First strategy: Add floor (i/2) and ceil (i/2). */
            j = i/2;
#endif
#if 1
            /* Second strategy:
               If i is even, add i/2.
               Otherwise, the previous strategy would add two numbers.
               Instead, add i-j for the largest occurring j less than i. */
            if (i % 2 == 0)
               j = i / 2;
            else
               for (j = i-1; b [j][0] == 0; j--);
#endif
#if 0
            /* Third strategy: As the second one, but for odd i, add
               an even i-j for the largest possible j; at the latest, this
               happens for i-j=1. */
            if (i % 2 == 0)
               j = i / 2;
            else
               for (j = i-1; j % 2 == 0 || b [j][0] == 0; j--);
#endif
            b [i][1] = j;
            if (b [j][0] == 0) {
               b [j][0] = 2;
               added++;
            }
            if (b [i-j][0] == 0) {
               b [i-j][0] = 2;
               added++;
            }
         }
      }

   return (added);
}

/*****************************************************************************/

static int minimal_dense_addition_sequence (long int ** b, int m)
   /* The parameters are the same as for dense_addition_sequence,
      except that b has three columns: b [i][2] is set to the desired
      precision if the exponent i occurs.
      The function returns an addition sequence that is minimal in the sense
      that none of the added entries may be removed (which does not mean that
      it is optimal!).
      Also, if possible it replaces remaining general additions by doublings
      (the addition sequence computation already privileges doublings, but
      these may become possible only later using added terms).
      The precision is also tracked, when writing an element as a sum of two
      smaller ones, then the precision of the smaller elements may need to be
      increased to that of the target. */
{
   int added, new_added, changed, i, j;
   long int **bnew;

   added = dense_addition_sequence (b, m);
   /* We need to work on a copy, since the computation of addition sequences
      via side effects destroys the previously computed sequence. It may
      happen that when removing one of the additional terms, three new ones
      are added in a run! Then we need to simply throw away the new
      sequence. */
   bnew = (long int **) malloc ((m + 1) * sizeof (long int *));
      for (i = 0; i <= m; i++)
         bnew [i] = (long *) malloc (2 * sizeof (long));

   /* This additional loop does not seem to be needed in practice. */
   changed = 1;
   while (changed) {
      changed = 0;
      /* Check if any of the added entries may be removed again. */
      for (i = m; i >= 2; i--)
         if (b [i][0] == 2) {
            /* Work on a copy of b with this element cancelled. */
            for (j = 0; j <= m; j++) {
               bnew [j][0] = b [j][0];
               bnew [j][1] = 0;
            }
            bnew [i][0] = 0;
            new_added = dense_addition_sequence (bnew, m);
            if (new_added == 0) {
               /* Copy the new addition sequence back. */
               for (j = 0; j <= m; j++) {
                  b [j][0] = bnew [j][0];
                  b [j][1] = bnew [j][1];
               }
               added--;
               changed = 1;
            }
         }
   }

   /* Look for doublings. */
   for (i = 4; i <= m; i += 2)
      if (b [i][0] != 0 && b [i/2][0] != 0)
         b [i][1] = i/2;

   /* Track the precision. */
   for (i = m; i >= 2; i--)
      if (b [i][0] != 0) {
         j = b [i][1];
         if (b [i][2] > b [j][2])
            b [j][2] = b [i][2];
         if (b [i][2] > b [i-j][2])
            b [i-j][2] = b [i][2];
      }

   for (i = 0; i <= m; i++)
      free (bnew [i]);
   free (bnew);

   return (added);
}

/*****************************************************************************/

static void qdev_eval_bsgs (ctype rop, cm_qdev_t f, ctype q1,
   double delta, int N)
   /* Evaluate f in q1 using the optimised addition sequence from f.
      N is the last index used in the addition chain.
      delta is the number of bits gained with each power of q.
      rop and q1 may be the same. */
{
   mp_prec_t prec, local_prec;
   int mopt [37] = { 2, 5, 7, 11, 13, 17, 19, 23, 55, 65, 77,
      91, 119, 133, 143, 175, 275, 325, 455, 595, 665, 715, 935, 1001,
      1309, 1463, 1547, 1729, 2275, 2975, 3325, 3575, 4675,
      6545, 7315, 7735, 8645 };
   int p [37] = { 2, 3, 4, 6, 7, 9, 10, 12, 18, 21, 24,
      28, 36, 40, 42, 44, 66, 77, 84, 108, 120, 126, 162, 168,
      216, 240, 252, 280, 308, 396, 440, 462, 594,
      648, 720, 756, 840 };
   long int T, cost, cost_new;
   int m, index;
   long int **bs;
   mpc_t *q, *c, tmp1, tmp2, tmp3;
   int i, j, k, J;

   prec = fget_prec (crealref (rop));
   cinit (tmp1, prec);
   cinit (tmp2, prec);
   cinit (tmp3, prec);
   T = f.chain [N][0];
   /* Find the optimal m=mopt[i] minimising the theoretical cost function
      (roughly, T/mopt[i] + p[i]). The function seems to be unimodular,
      so we take the smallest i before it increases again. */
   i = 0;
   cost_new = (T + mopt [0]) / mopt[0] + p[0] - 1;
   do {
      cost = cost_new;
      i++;
      cost_new = (T + mopt[i]) / mopt[i] + p[i] - 1;
   } while (cost_new < cost && i < 36);
   if (i == 36 && cost_new < cost) {
      printf ("*** Houston, we have a problem!\n");
      printf ("mopt and p too short in 'qdev_eval_bsgs'.\n");
      exit (1);
   }
   else
      m = mopt [i-1];

   /* Determine the occurring baby-steps; in practice, these are as
      many as predicted by p. Also keep track of their required precision. */
   bs = (long int **) malloc ((m + 1) * sizeof (long int *));
   for (i = 0; i <= m; i++) {
      bs [i] = (long *) malloc (3 * sizeof (long));
      for (j = 0; j < 3; j++)
         bs [i][j] = 0;
   }
   for (i = 0; i <= N; i++) {
      index = f.chain [i][0] % m;
      /* Register the precision needed for the term with lowest exponent
         and not later ones. */
      if (bs [index][0] == 0) {
         bs [index][0] = 1;
         bs [index][2] = (long int) prec
                         - (long int) (f.chain [i][0] * delta);
      }
   }
   bs [m][0] = 1; /* for the giant steps */
   bs [m][2] = (long int) prec - (long int) (m * delta);
   minimal_dense_addition_sequence (bs, m);

   /* Compute the baby-steps. */
   q = (mpc_t *) malloc ((m + 1) * sizeof (mpc_t));
   cinit (q [0], 2);
   cset_ui (q [0], 1);
   cinit (q [1], bs [1][2]);
   cset (q [1], q1);
   for (i = 2; i <= m; i++)
      if (bs [i][0] != 0) {
         cinit (q [i], bs [i][2]);
         j = bs [i][1];
         k = i - j;
         /* Here decreasing the argument precision to the target precision
            apparently does not gain time; both are probably too close for
            it to make a difference. */
         if (j == k)
            csqr (q [i], q [j]);
         else
            cmul (q [i], q [j], q [k]);
      }

   /* Compute the giant steps; we need
      \sum_{j=0}^{J-1} (\sum_{k=0}^{m-1} c_{k+j*m} q^k) (q^m)^j.
      First compute the inner coefficients, then use a Horner scheme. */
   J = f.chain [N][0] / m + 1;
   c = (mpc_t *) malloc (J * sizeof (mpc_t));
   for (j = 0; j < J; j++) {
      local_prec = prec - (mp_prec_t) (j * m * delta);
      cinit (c [j], local_prec);
      cset_ui (c [j], 0);
   }
   for (i = 0; i <= N; i++)
      if (f.chain [i][4] != 0) {
         j = f.chain [i][0] / m;
         k = f.chain [i][0] % m; 
         /* We assume the coefficients are 1 or -1. */
         if (f.chain [i][4] == 1)
            cadd (c [j], c[j], q [k]);
         else
            csub (c [j], c[j], q [k]);
      }
   cset (rop, c [J-1]);
   for (j = J-2; j >= 0; j--) {
      /* Carry out the multiplication at the precision of c [j+1]. */
      local_prec = prec - (mp_prec_t) ((j+1) * m * delta);
      cset_prec (tmp1, local_prec);
      cset_prec (tmp2, local_prec);
      cset_prec (tmp3, local_prec);
      cset (tmp1, rop);
      cset (tmp2, q [m]);
      cmul (tmp3, tmp1, tmp2);
      cadd (rop, tmp3, c [j]);
   }

   for (i = 0; i <= m; i++) {
      if (bs [i][0] != 0)
         cclear (q [i]);
      free (bs [i]);
   }
   for (j = 0; j < J; j++)
      cclear (c [j]);
   free (bs);
   free (q);
   free (c);
   cclear (tmp1);
   cclear (tmp2);
   cclear (tmp3);
}

/*****************************************************************************/

void cm_qdev_eval (ctype rop, cm_qdev_t f, ctype q1)
   /* Evaluate f in q1. rop and q1 may be the same. */

{
   mp_prec_t prec;
   double    delta;
   long int T;
   int N;

   /* Compute the last exponent T and its index N in f.chain. */
   prec = fget_prec (crealref (rop));
   delta = - lognorm2 (q1);
   T = (prec - 2) / delta;
   for (N=0; N < f.length && f.chain [N][0] <= T; N++);
   if (N == f.length) {
      printf ("*** Houston, we have a problem! Addition chain too short ");
      printf ("in 'cm_qdev_eval'.\n");
      printf ("T=%li, length=%i\n", T, f.length);
      exit (1);
   }
   N--;
   T = f.chain [N][0];

   if (N < 20)
      qdev_eval_addition_sequence (rop, f, q1, delta, N);
   else
      qdev_eval_bsgs (rop, f, q1, delta, N);
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

