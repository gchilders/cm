
/*

mpfpx.c - mpfpx library

Copyright (C) 2009 Andreas Enge

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

#include "mpfpx-impl.h"

#define MAX(a,b) (((a) > (b) ? (a) : (b)))
#define MIN(a,b) (((a) < (b) ? (a) : (b)))

static int __mpfpx_arithmetic;

/*****************************************************************************/

static void mpfpx_karatsuba_split (mpfpx_t a, mpfpx_t b, mpfpx_t op)
   /* writes op as a (X^2) X + b (X^2), where a, b and op must be distinct */

{
   int i;

   a->deg = (op->deg - 1) / 2;
   b->deg = op->deg / 2;
   for (i = 0; i <= a->deg; i++)
   {
      mpfp_set (a->coeff [i], op->coeff [2*i+1]);
      mpfp_set (b->coeff [i], op->coeff [2*i]);
   }
   if (a->deg != b->deg)
      mpfp_set (b->coeff [b->deg], op->coeff [2*b->deg]);
}

/*****************************************************************************/

static void mpfpx_newton (mpfpx_t rop, mpfpx_t op, int exp)
   /* computes rop such that rop*op \equiv 1 \mod X^exp; see Algorithm 9.3 */
   /* in GG99                                                              */

{
   int prec = 1;
      /* The result is exact modulo X^prec. */
   mpfpx_t lop, tmp;

   mpfpx_init_set (lop, op);
   mpfpx_init (tmp);

   rop->deg = 0;
   mpfp_inv (rop->coeff [0], lop->coeff [0]);

   while (prec < exp)
   {
      prec = MIN(2 * prec, exp);
      mpfpx_sqr (tmp, rop);
      mpfpx_mul (tmp, tmp, lop);
      mpfpx_sub (tmp, rop, tmp);
      mpfpx_add (rop, rop, tmp);
      rop->deg = prec - 1;
   }

   mpfpx_clear (lop);
   mpfpx_clear (tmp);
}

/*****************************************************************************/

static void mpfpx_rev (mpfpx_t rop, mpfpx_t op, int k)
   /* puts rop = X^k * op (1/X); only meaningful for k >= deg op */

{
   int     i;
   mpfpx_t lrop;

   mpfpx_init (lrop);

   for (i = 0; i <= op->deg; i++)
      mpfp_set (lrop->coeff [k - i], op->coeff [i]);
   for (i = 0; i < k - op->deg; i++)
      mpfp_set_ui (lrop->coeff [i], 0ul);
   lrop->deg = k;
   mpfpx_normalise (lrop);
   mpfpx_set (rop, lrop);

   mpfpx_clear (lrop);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_type_init (mpz_t p, int arithmetic)

{
   mpfp_type_init (p);
   __mpfpx_arithmetic = arithmetic;
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_init (mpfpx_t rop)

{
   int i;

   for (i = 0; i <= MPFPX_MAX_DEG; i++)
      mpfp_init (rop->coeff [i]);
   mpfp_set_ui (rop->coeff [0], 0ul);
   rop->deg = -1;
}

/*****************************************************************************/

void mpfpx_clear (mpfpx_t rop)

{
   int i;

   for (i = 0; i <= MPFPX_MAX_DEG; i++)
      mpfp_clear (rop->coeff [i]);
}

/*****************************************************************************/

void mpfpx_normalise (mpfpx_t rop)
   /* deletes leading zeroes */

{
   while (rop->deg > 0 && mpfp_cmp_ui (rop->coeff [rop->deg], 0ul) == 0)
      rop->deg --;
   if (rop->deg == 0 && mpfp_cmp_ui (rop->coeff [rop->deg], 0ul) == 0)
      rop->deg = -1;
}

/*****************************************************************************/

void mpfpx_monicise (mpfpx_t rop, mpfpx_t op, mpfp_t lc_inv)

{
   int i;

   rop->deg = op->deg;
   if (lc_inv == NULL)
   {
      mpfp_t inv;
      mpfp_init (inv);
      mpfp_inv (inv, op->coeff [op->deg]);
      for (i = 0; i < rop->deg; i++)
         mpfp_mul (rop->coeff [i], op->coeff [i], inv);
      mpfp_set_ui (rop->coeff [op->deg], 1ul);
      mpfp_clear (inv);
   }
   else
   {
      for (i = 0; i < rop->deg; i++)
         mpfp_mul (rop->coeff [i], op->coeff [i], lc_inv);
      mpfp_set_ui (rop->coeff [op->deg], 1ul);
   }
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_set (mpfpx_t rop, mpfpx_t op)

{
   int i;

   rop->deg = op->deg;

   for (i = 0; i <= op->deg; i++)
      mpfp_set (rop->coeff [i], op->coeff [i]);
}

/*****************************************************************************/

void mpfpx_set_fp (mpfpx_t rop, mpfp_t op)

{
   if (mpfp_cmp_ui (op, 0ul) == 0)
      rop->deg = -1;
   else
   {
      rop->deg = 0;
      mpfp_set (rop->coeff [0], op);
   }
}

/*****************************************************************************/

void mpfpx_set_ui (mpfpx_t rop, unsigned long int op)

{
   if (op == 0)
      rop->deg = -1;
   else
   {
      rop->deg = 0;
      mpfp_set_ui (rop->coeff [0], op);
   }
}

/*****************************************************************************/

void mpfpx_set_z_array (mpfpx_t rop, mpz_t *op, int size)

{
   int i;

   rop->deg = size - 1;
   for (i = 0; i <= rop->deg; i++)
      mpfp_set_z (rop->coeff [i], op [i]);
}

/*****************************************************************************/

void mpfpx_set_ui_array (mpfpx_t rop, const unsigned long int *op, int size)

{
   int i;

   rop->deg = size - 1;
   for (i = 0; i <= rop->deg; i++)
      mpfp_set_ui (rop->coeff [i], op [i]);
}

/*****************************************************************************/

void mpfpx_set_str_array (mpfpx_t rop, char **op, int base, int size)

{
   int i;

   rop->deg = size - 1;
   for (i = 0; i <= rop->deg; i++)
      mpfp_set_str (rop->coeff [i], op [i], base);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_init_set (mpfpx_t rop, mpfpx_t op)

{
   mpfpx_init (rop);
   mpfpx_set (rop, op);
}

/*****************************************************************************/

void mpfpx_init_set_fp (mpfpx_t rop, mpfp_t op)

{
   mpfpx_init (rop);
   mpfpx_set_fp (rop, op);
}

/*****************************************************************************/

void mpfpx_init_set_ui (mpfpx_t rop, unsigned long op)

{
   mpfpx_init (rop);
   mpfpx_set_ui (rop, op);
}

/*****************************************************************************/

void mpfpx_init_set_ui_array (mpfpx_t rop, const unsigned long int *op,
   int size)

{
   mpfpx_init (rop);
   mpfpx_set_ui_array (rop, op, size);
}

/*****************************************************************************/

void mpfpx_init_set_str_array (mpfpx_t rop, char **op, int base, int size)

{
   mpfpx_init (rop);
   mpfpx_set_str_array (rop, op, base, size);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_out (mpfpx_t op)

{
   int i;

   printf ("[");
   if (op->deg == -1)
      printf ("*");
   else
   {
      for (i = 0; i < op->deg; i++)
      {
         mpfp_out (op->coeff [i]);
         printf (" ");
      }
      mpfp_out (op->coeff [op->deg]);
   }
   printf ("]\n");
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_deg_set (mpfpx_t rop, int op)

{
   rop->deg = op;
}

/*****************************************************************************/

void mpfpx_coeff_set_z (mpfpx_t rop, unsigned int op1, mpz_t op2)

{
   mpfp_set_z (rop->coeff [op1], op2);
}

/*****************************************************************************/

void mpfpx_coeff_set_ui (mpfpx_t rop, unsigned int op1, unsigned long int op2)

{
   mpfp_set_ui (rop->coeff [op1], op2);
}

/*****************************************************************************/

void mpfpx_coeff_get_z (mpz_t rop, mpfpx_t op1, unsigned int op2)

{
   mpfp_get_z (rop, op1->coeff [op2]);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_add (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2)

{
   __mpfpx_t *lop1, *lop2;
   int        i;

   if (op1->deg >= op2->deg)
   {
      lop1 = op1;
      lop2 = op2;
   }
   else
   {
      lop1 = op2;
      lop2 = op1;
   }

   for (i = 0; i <= lop2->deg; i++)
      mpfp_add (rop->coeff [i], lop1->coeff [i], lop2->coeff [i]);
   for (i = lop2->deg + 1; i <= lop1->deg; i++)
      mpfp_set (rop->coeff [i], lop1->coeff [i]);
   rop->deg = lop1->deg;
}

/*****************************************************************************/

void mpfpx_sub (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2)

{
   mpfpx_t tmp;

   mpfpx_init (tmp);

   mpfpx_neg (tmp, op2);
   mpfpx_add (rop, op1, tmp);

   mpfpx_clear (tmp);
}

/*****************************************************************************/

void mpfpx_sub_ui (mpfpx_t rop, mpfpx_t op1, unsigned long int op2)

{
   int i;

   if (op1->deg == -1)
      if (op2 == 0)
         rop->deg = -1;
      else
      {
         rop->deg = 0;
         mpfp_set_ui (rop->coeff [0], op2);
         mpfp_neg (rop->coeff [0], rop->coeff [0]);
      }
   else if (rop == op1)
      mpfp_sub_ui (rop->coeff [0], rop->coeff [0], op2);
   else
   {
      rop->deg = op1->deg;
      mpfp_sub_ui (rop->coeff [0], rop->coeff [0], op2);
      for (i = 1; i <= op1->deg; i++)
         mpfp_set (rop->coeff [i], op1->coeff [i]);
   }
}

/*****************************************************************************/

void mpfpx_mul (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2)

{
   int     i;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op1->deg == 0)
   {
      if (mpfp_cmp_ui (op1->coeff [0], 1ul) == 0)
         mpfpx_set (rop, op2);
      else if (mpfp_cmp_si (op1->coeff [0], -1l) == 0)
         mpfpx_neg (rop, op2);
      else
      {
         rop->deg = op2->deg;
         for (i = 0; i <= op2->deg; i++)
            mpfp_mul (rop->coeff [i], op2->coeff [i], op1->coeff [0]);
      }
   }
   else if (op2->deg == 0)
   {
      if (mpfp_cmp_ui (op2->coeff [0], 1ul) == 0)
         mpfpx_set (rop, op1);
      else if (mpfp_cmp_si (op2->coeff [0], -1l) == 0)
         mpfpx_neg (rop, op1);
      else
      {
         rop->deg = op1->deg;
         for (i = 0; i <= op1->deg; i++)
            mpfp_mul (rop->coeff [i], op1->coeff [i], op2->coeff [0]);
      }
   }
   else
   if (__mpfpx_arithmetic == MPFPX_KARATSUBA)
   {
      /* We put op1 = a (X^2) X + b (X^2) and op2 = c (X^2) X + d (X^2) */
      /* and do plenty of copying                                       */
      mpfpx_t a, b, c, d, u, v, w;

      mpfpx_init (a);
      mpfpx_init (b);
      mpfpx_init (c);
      mpfpx_init (d);
      mpfpx_init (u);
      mpfpx_init (v);
      mpfpx_init (w);

      mpfpx_karatsuba_split (a, b, op1);
      mpfpx_karatsuba_split (c, d, op2);
      mpfpx_add (v, a, b);
      mpfpx_add (w, c, d);
      mpfpx_mul (u, v, w);
      mpfpx_mul (v, a, c);
      mpfpx_mul (w, b, d);
      rop->deg = op1->deg + op2->deg;
      /* we have to put rop = v (X^2) X^2 +  w (X^2) + (u - v - w)(X^2) X */
      mpfpx_add (a, v, w);
      mpfpx_sub (b, u, a);
      /* copy b (X^2) X and w (X^2) and clear further coefficients */
      for (i = 0; i <= b->deg; i++)
         mpfp_set (rop->coeff [2*i+1], b->coeff [i]);
      for (i = 2 * b->deg + 3; i <= rop->deg; i += 2);
         mpfp_set_ui (rop->coeff [i], 0ul);
      for (i = 0; i <= w->deg; i++)
         mpfp_set (rop->coeff [2*i], w->coeff [i]);
      for (i = 2 * w->deg + 2; i <= rop->deg; i += 2)
         mpfp_set_ui (rop->coeff [i], 0ul);
      /* add v (X^2) X^2 */
      for (i = 0; i <= v->deg; i++)
         mpfp_add (rop->coeff [2*i+2], rop->coeff [2*i+2], v->coeff [i]);

      mpfpx_clear (a);
      mpfpx_clear (b);
      mpfpx_clear (c);
      mpfpx_clear (d);
      mpfpx_clear (u);
      mpfpx_clear (v);
      mpfpx_clear (w);
/*
      // We put op1 = a (X) X^n + b (X) and op2 = c (X) X^n + d (X)
      // and do plenty of copying
      mpfpx_t a, b, c, d, u, v, w;
      int     n = max (op1-> deg / 2, op2->deg / 2) + 1;

      mpfpx_init (a);
      mpfpx_init (b);
      mpfpx_init (c);
      mpfpx_init (d);
      mpfpx_init (u);
      mpfpx_init (v);
      mpfpx_init (w);

      mpfpx_karatsuba_split2 (a, b, op1, n);
      mpfpx_karatsuba_split2 (c, d, op2, n);
      mpfpx_add (v, a, b);
      mpfpx_add (w, c, d);
      mpfpx_mul (u, v, w);
      mpfpx_mul (v, a, c);
      mpfpx_mul (w, b, d);
      rop->deg = op1->deg + op2->deg;
      // we have to put rop = v (X) X^{2n} + (u - v - w)(X) X^n + w (X)
      // copy v (X) X^{2n} and w (X) and clear further coefficients
      for (i = 0; i <= w->deg; i++)
         mpfp_set (rop->coeff [i], w->coeff [i]);
      for (i = w->deg + 1; i < 2*n; i++)
         mpfp_set_ui (rop->coeff [i], 0);
      for (i = 0; i <= v->deg; i++)
         mpfp_set (rop->coeff [i + 2*n], v->coeff [i]);
      for (i = 2*n + v->deg + 1; i <= rop->deg; i++)
         mpfp_set_ui (rop->coeff [i], 0);
      mpfpx_add (a, v, w);
      mpfpx_sub (b, u, a);
      // add b (X) X^n
      for (i = 0; i <= b->deg; i++)
         mpfp_add (rop->coeff [i+n], rop->coeff [i+n], b->coeff [i]);
      mpfpx_normalise (rop);

      mpfpx_clear (a);
      mpfpx_clear (b);
      mpfpx_clear (c);
      mpfpx_clear (d);
      mpfpx_clear (u);
      mpfpx_clear (v);
      mpfpx_clear (w);
*/
   }
   else
   {
      /* naive arithmetic */
      mpfpx_t lrop;
      int     j;

      mpfpx_init (lrop);

      lrop->deg = op1->deg + op2->deg;
      for (i = 0; i <= lrop->deg; i++)
         mpfp_set_ui (lrop->coeff [i], 0ul);

      for (i = 0; i <= op1->deg; i++)
         for (j = 0; j <= op2->deg; j++)
            mpfp_add_mul (lrop->coeff [i+j], op1->coeff [i], op2->coeff [j]);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }

}

/*****************************************************************************/

void mpfpx_mul_ui (mpfpx_t rop, mpfpx_t op1, unsigned long int op2)

{
   int     i;

   rop->deg = op1->deg;
   for (i = 0; i <= op1->deg; i++)
       mpfp_mul_ui (rop->coeff [i], op1->coeff [i], op2);
}

/*****************************************************************************/

void mpfpx_mul_low (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp)

{
   int     i;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op1->deg == 0)
   {
      if (mpfp_cmp_ui (op1->coeff [0], 1ul) == 0)
         mpfpx_set (rop, op2);
      else if (mpfp_cmp_si (op1->coeff [0], -1l) == 0)
         mpfpx_neg (rop, op2);
      else
      {
         rop->deg = op2->deg;
         for (i = 0; i <= MIN(op2->deg, exp); i++)
            mpfp_mul (rop->coeff [i], op2->coeff [i], op1->coeff [0]);
      }
   }
   else if (op2->deg == 0)
   {
      if (mpfp_cmp_ui (op2->coeff [0], 1ul) == 0)
         mpfpx_set (rop, op1);
      else if (mpfp_cmp_si (op2->coeff [0], -1l) == 0)
         mpfpx_neg (rop, op1);
      else
      {
         rop->deg = op1->deg;
         for (i = 0; i <= MIN(op1->deg, exp); i++)
            mpfp_mul (rop->coeff [i], op1->coeff [i], op2->coeff [0]);
      }
   }
   else
   if (__mpfpx_arithmetic == MPFPX_KARATSUBA && exp >= 3)
   {
      int op1_deg = op1->deg, op2_deg = op2->deg;

      op1->deg = MIN(op1->deg, exp);
      op2->deg = MIN(op2->deg, exp);
      mpfpx_mul (rop, op1, op2);
      /* resetting degrees in the following order in case rop is one of op1 */
      /* or op2                                                             */
      op1->deg = op1_deg;
      op2->deg = op2_deg;
      rop->deg = op1_deg + op2_deg;
   }
   else
   {
      /* naive arithmetic */
      mpfpx_t lrop;
      int     j;

      mpfpx_init (lrop);

      lrop->deg = op1->deg + op2->deg;
      for (i = 0; i <= lrop->deg; i++)
         mpfp_set_ui (lrop->coeff [i], 0ul);

      for (i = 0; i <= MIN(op1->deg, exp); i++)
         for (j = 0; j <= MIN(op2->deg, exp - i); j++)
            mpfp_add_mul (lrop->coeff [i+j], op1->coeff [i], op2->coeff [j]);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }

}

/*****************************************************************************/

void mpfpx_mul_high (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp)

{
   int i, j;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op1->deg == 0) {
      if (mpfp_cmp_ui (op1->coeff [0], 1ul) == 0)
         mpfpx_set (rop, op2);
      else if (mpfp_cmp_si (op1->coeff [0], -1l) == 0)
         mpfpx_neg (rop, op2);
      else {
         rop->deg = op2->deg;
         for (i = op2->deg; i >= MAX(op2->deg - exp, 0); i--)
            mpfp_mul (rop->coeff [i], op2->coeff [i], op1->coeff [0]);
      }
   }
   else if (op2->deg == 0) {
      if (mpfp_cmp_ui (op2->coeff [0], 1ul) == 0)
         mpfpx_set (rop, op1);
      else if (mpfp_cmp_si (op2->coeff [0], -1l) == 0)
         mpfpx_neg (rop, op1);
      else {
         rop->deg = op1->deg;
         for (i = op1->deg; i >= MAX(op1->deg - exp, 0); i--)
            mpfp_mul (rop->coeff [i], op1->coeff [i], op2->coeff [0]);
      }
   }
   else if (__mpfpx_arithmetic == MPFPX_KARATSUBA && exp >= 3) {
      /* copying the first coefficients of the operators and then using */
      /* Karatsuba, which is probably not optimal                       */
      mpfpx_t lop1, lop2;
      int     op1_deg, op2_deg;

      mpfpx_init (lop1);
      mpfpx_init (lop2);

      op1_deg = op1->deg;
      op2_deg = op2->deg;
      lop1->deg = MIN(op1_deg, exp);
      for (i = 0; i <= lop1->deg; i++)
         mpfp_set (lop1->coeff [lop1->deg - i], op1->coeff [op1->deg - i]);
      lop2->deg = MIN(op2_deg, exp);
      for (i = 0; i <= lop2->deg; i++)
         mpfp_set (lop2->coeff [lop2->deg - i], op2->coeff [op2->deg - i]);
      mpfpx_mul (rop, lop1, lop2);
      for (i = 0; i <= MIN(exp, rop->deg); i++)
         mpfp_set (rop->coeff [op1_deg + op2_deg - i], rop->coeff [rop->deg - i]);
      rop->deg = op1_deg + op2_deg;

      mpfpx_clear (lop1);
      mpfpx_clear (lop2);
   }
   else {
      /* naive arithmetic */
      mpfpx_t lrop;

      mpfpx_init (lrop);

      lrop->deg = op1->deg + op2->deg;
      for (i = 0; i <= lrop->deg; i++)
         mpfp_set_ui (lrop->coeff [i], 0ul);

      for (i = op1->deg; i >= MAX(op1->deg - exp, 0); i--)
         for (j = op2->deg; j >= MAX(lrop->deg - exp - i, 0); j--)
            mpfp_add_mul (lrop->coeff [i+j], op1->coeff [i], op2->coeff [j]);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }
}

/*****************************************************************************/

void mpfpx_mul_m (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2)

{
   mpfpx_t lrop;
   int     i;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op2->deg == 0)
      mpfpx_set (rop, op1);
   else
   {
      mpfpx_init (lrop);
      op2->deg --;
      mpfpx_mul (lrop, op1, op2);
      op2->deg ++;
      lrop->deg = op1->deg + op2->deg;
      for (i = 0; i <= op1->deg; i++)
         mpfp_add (lrop->coeff [i + op2->deg], lrop->coeff [i + op2->deg],
                   op1->coeff [i]);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }
}

/*****************************************************************************/

void mpfpx_mul_m_low (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp)

{
   mpfpx_t lrop;
   int     i;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op2->deg == 0)
      mpfpx_set (rop, op1);
   else
   {
      mpfpx_init (lrop);
      op2->deg --;
      mpfpx_mul_low (lrop, op1, op2, exp);
      op2->deg ++;
      lrop->deg = op1->deg + op2->deg;
      for (i = 0; i <= MIN(op1->deg, exp); i++)
         mpfp_add (lrop->coeff [i + op2->deg], lrop->coeff [i + op2->deg],
                   op1->coeff [i]);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }
}

/*****************************************************************************/

void mpfpx_mul_m_high (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp)

{
   mpfpx_t lrop;
   int     i;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op2->deg == 0)
      mpfpx_set (rop, op1);
   else
   {
      mpfpx_init (lrop);
      op2->deg --;
      mpfpx_mul_high (lrop, op1, op2, exp - 1);
      op2->deg ++;
      lrop->deg = op1->deg + op2->deg;
      for (i = MAX(op1->deg - exp, 0); i <= op1->deg; i++)
         mpfp_add (lrop->coeff [i + op2->deg], lrop->coeff [i + op2->deg],
                   op1->coeff [i]);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }
}

/*****************************************************************************/

void mpfpx_mul_mm_high (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp)

{
   mpfpx_t lrop;
   int     i;

   if (op1->deg == -1 || op2->deg == -1)
      rop->deg = -1;
   else if (op2->deg == 0)
      mpfpx_set (rop, op1);
   else if (op1->deg == 0)
      mpfpx_set (rop, op2);
   else
   {
      mpfpx_init (lrop);
      op1->deg --;
      op2->deg --;
      mpfpx_mul_high (lrop, op1, op2, exp - 2);
      op1->deg ++;
      op2->deg ++;
      lrop->deg = op1->deg + op2->deg;
      for (i = MAX(op1->deg - exp, 0); i < op1->deg; i++)
         mpfp_add (lrop->coeff [i + op2->deg], lrop->coeff [i + op2->deg],
                   op1->coeff [i]);
      for (i = MAX(op2->deg - exp, 0); i < op2->deg; i++)
         mpfp_add (lrop->coeff [i + op1->deg], lrop->coeff [i + op1->deg],
                   op2->coeff [i]);
      mpfp_set_ui (lrop->coeff [lrop->deg], 1ul);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }
}

/*****************************************************************************/

void mpfpx_sqr (mpfpx_t rop, mpfpx_t op)

{
   int     i, j;

   if (op->deg == -1)
      rop->deg = -1;
   else if (op->deg == 0)
   {
      rop->deg = 0;
      if (   mpfp_cmp_ui (op->coeff [0],  1ul) == 0
          || mpfp_cmp_si (op->coeff [0], -1l) == 0)
         mpfp_set_ui (rop->coeff [0], 1ul);
      else
         mpfp_sqr (rop->coeff [0], op->coeff [0]);
   }
   else
   if (__mpfpx_arithmetic == MPFPX_KARATSUBA)
   {
      if (op->deg == 2)
      {
         /* We need a special treatment of the case of degree 2, where          */
         /* Karatsuba's approach would need 7 squarings, while naive            */
         /* multiplication needs 3 squarings and 3 multiplications. We can make */
         /* do with 6 squarings.                                                */
         mpfp_t u, v, w, x, y, z;

         mpfp_init (u);
         mpfp_init (v);
         mpfp_init (w);
         mpfp_init (x);
         mpfp_init (y);
         mpfp_init (z);

         rop->deg = 4;
         mpfp_add (x, op->coeff [2], op->coeff [1]);
         mpfp_sqr (u, x);
         mpfp_add (x, op->coeff [2], op->coeff [0]);
         mpfp_sqr (v, x);
         mpfp_add (x, op->coeff [1], op->coeff [0]);
         mpfp_sqr (w, x);
         mpfp_sqr (x, op->coeff [2]);
         mpfp_sqr (y, op->coeff [1]);
         mpfp_sqr (z, op->coeff [0]);
         mpfp_set (rop->coeff [4], x);
         mpfp_set (rop->coeff [3], u);
         mpfp_sub (rop->coeff [3], rop->coeff [3], x);
         mpfp_sub (rop->coeff [3], rop->coeff [3], y);
         mpfp_set (rop->coeff [2], v);
         mpfp_sub (rop->coeff [2], rop->coeff [2], x);
         mpfp_sub (rop->coeff [2], rop->coeff [2], z);
         mpfp_add (rop->coeff [2], rop->coeff [2], y);
         mpfp_set (rop->coeff [1], w);
         mpfp_sub (rop->coeff [1], rop->coeff [1], y);
         mpfp_sub (rop->coeff [1], rop->coeff [1], z);
         mpfp_set (rop->coeff [0], z);

         mpfp_clear (u);
         mpfp_clear (v);
         mpfp_clear (w);
         mpfp_clear (x);
         mpfp_clear (y);
         mpfp_clear (z);

      }
      else
      {
         /* for details, see mpfpx_mul */
         mpfpx_t a, b, u, v, w;

         mpfpx_init (a);
         mpfpx_init (b);
         mpfpx_init (u);
         mpfpx_init (v);
         mpfpx_init (w);

         mpfpx_karatsuba_split (a, b, op);
         mpfpx_add (v, a, b);
         mpfpx_sqr (u, v);
         mpfpx_sqr (v, a);
         mpfpx_sqr (w, b);
         rop->deg = 2 * op->deg;
         mpfpx_add (a, v, w);
         mpfpx_sub (b, u, a);
         for (i = 0; i <= b->deg; i++)
            mpfp_set (rop->coeff [2*i+1], b->coeff [i]);
         for (i = 2 * b->deg + 3; i <= rop->deg; i += 2);
            mpfp_set_ui (rop->coeff [i], 0ul);
         for (i = 0; i <= w->deg; i++)
            mpfp_set (rop->coeff [2*i], w->coeff [i]);
         for (i = 2 * w->deg + 2; i <= rop->deg; i += 2)
            mpfp_set_ui (rop->coeff [i], 0ul);
         for (i = 0; i <= v->deg; i++)
            mpfp_add (rop->coeff [2*i+2], rop->coeff [2*i+2], v->coeff [i]);

         mpfpx_clear (a);
         mpfpx_clear (b);
         mpfpx_clear (u);
         mpfpx_clear (v);
         mpfpx_clear (w);
      }
   }
   else
   {
      /* naive arithmetic */
      mpfpx_t lrop;
      mpfp_t  tmp;

      mpfpx_init (lrop);
      mpfp_init (tmp);

      lrop->deg = 2 * op->deg;
      for (i = 0; i <= lrop->deg; i++)
         mpfp_set_ui (lrop->coeff [i], 0ul);

      for (i = 0; i <= op->deg; i++)
      {
         mpfp_sqr (tmp, op->coeff [i]);
         mpfp_add (lrop->coeff [2*i], lrop->coeff [2*i], tmp);
         for (j = i+1; j <= op->deg; j++)
         {
            mpfp_mul (tmp, op->coeff [i], op->coeff [j]);
            mpfp_add (lrop->coeff [i+j], lrop->coeff [i+j], tmp);
            mpfp_add (lrop->coeff [i+j], lrop->coeff [i+j], tmp);
         }
      }

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
      mpfp_clear (tmp);
   }
}

/*****************************************************************************/

void mpfpx_sqr_low (mpfpx_t rop, mpfpx_t op, char exp)

{
   int     i, j;

   if (op->deg == -1)
      rop->deg = -1;
   else if (op->deg == 0)
   {
      rop->deg = 0;
      if (   mpfp_cmp_ui (op->coeff [0],  1ul) == 0
          || mpfp_cmp_si (op->coeff [0], -1l) == 0)
         mpfp_set_ui (rop->coeff [0], 1ul);
      else
         mpfp_sqr (rop->coeff [0], op->coeff [0]);
   }
   else
   if (__mpfpx_arithmetic == MPFPX_KARATSUBA)
   {
      printf ("\n*** Trying Karatsuba in mpfpx_sqr_low!\n");
      exit (1);
   }
   else
   {
      /* naive arithmetic */
      mpfpx_t lrop;
      mpfp_t  tmp;

      mpfpx_init (lrop);
      mpfp_init (tmp);

      lrop->deg = 2 * op->deg;
      for (i = 0; i <= lrop->deg; i++)
         mpfp_set_ui (lrop->coeff [i], 0ul);

      for (i = 0; i <= MIN(op->deg, exp / 2); i++)
      {
         mpfp_sqr (tmp, op->coeff [i]);
         mpfp_add (lrop->coeff [2*i], lrop->coeff [2*i], tmp);
      }
      for (i = 0; i <= MIN(op->deg, exp); i++)
         for (j = i+1; j <= MIN(op->deg, exp - i); j++)
         {
            mpfp_mul (tmp, op->coeff [i], op->coeff [j]);
            mpfp_add (lrop->coeff [i+j], lrop->coeff [i+j], tmp);
            mpfp_add (lrop->coeff [i+j], lrop->coeff [i+j], tmp);
         }

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
      mpfp_clear (tmp);
   }
}

/*****************************************************************************/

void mpfpx_sqr_high (mpfpx_t rop, mpfpx_t op, char exp)

{
   int     i, j;

   if (op->deg == -1)
      rop->deg = -1;
   else if (op->deg == 0) {
      rop->deg = 0;
      if (   mpfp_cmp_ui (op->coeff [0],  1ul) == 0
          || mpfp_cmp_si (op->coeff [0], -1l) == 0)
         mpfp_set_ui (rop->coeff [0], 1ul);
      else
         mpfp_sqr (rop->coeff [0], op->coeff [0]);
   }
   else
   if (__mpfpx_arithmetic == MPFPX_KARATSUBA && exp >= 5) {
      /* copying the first coefficients of the operator and then using       */
      /* Karatsuba, which is probably not optimal. The degree bound of 5 has */
      /* been determined experimentally; for 4, we need more multiplications */
      /* altogether (6 additional squares versus 5 saved multiplications)    */
      mpfpx_t lop;
      int     op_deg;

      mpfpx_init (lop);

      op_deg = op->deg;
      lop->deg = MIN(op_deg, exp);
      for (i = 0; i <= lop->deg; i++)
         mpfp_set (lop->coeff [lop->deg - i], op->coeff [op->deg - i]);
      mpfpx_sqr (rop, lop);
      for (i = 0; i <= MIN(exp, rop->deg); i++)
         mpfp_set (rop->coeff [2 * op_deg - i], rop->coeff [rop->deg - i]);
      rop->deg = 2 * op_deg;

      mpfpx_clear (lop);
   }
   else {
      /* naive arithmetic */
      mpfpx_t lrop;
      mpfp_t  tmp;

      mpfpx_init (lrop);
      mpfp_init (tmp);

      lrop->deg = 2 * op->deg;
      for (i = 0; i <= lrop->deg; i++)
         mpfp_set_ui (lrop->coeff [i], 0ul);

      for (i = op->deg; i >= MAX((lrop->deg - exp + 1) / 2, 0); i--) {
         mpfp_sqr (tmp, op->coeff [i]);
         mpfp_add (lrop->coeff [2*i], lrop->coeff [2*i], tmp);
      }
      for (i = op->deg; i >= MAX((lrop->deg - exp + 2) / 2, 0); i--) {
         for (j = i-1; j >= MAX(lrop->deg - exp - i, 0); j--) {
            mpfp_mul (tmp, op->coeff [i], op->coeff [j]);
            mpfp_add (lrop->coeff [i+j], lrop->coeff [i+j], tmp);
            mpfp_add (lrop->coeff [i+j], lrop->coeff [i+j], tmp);
         }
      }

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
      mpfp_clear (tmp);
   }
}

/*****************************************************************************/

void mpfpx_sqr_m_high (mpfpx_t rop, mpfpx_t op, char exp)

{
   mpfpx_t lrop;
   int     i;

   if (op->deg == -1)
      rop->deg = -1;
   else if (op->deg == 0)
      mpfpx_set (rop, op);
   else
   {
      mpfpx_init (lrop);
      op->deg --;
      mpfpx_sqr_high (lrop, op, exp - 2);
      op->deg ++;
      lrop->deg = 2 * op->deg;
      for (i = MAX(op->deg - exp, 0); i < op->deg; i++)
      {
         mpfp_add (lrop->coeff [i + op->deg], lrop->coeff [i + op->deg],
                   op->coeff [i]);
         mpfp_add (lrop->coeff [i + op->deg], lrop->coeff [i + op->deg],
                   op->coeff [i]);
      }
      mpfp_set_ui (lrop->coeff [lrop->deg], 1ul);

      mpfpx_set (rop, lrop);

      mpfpx_clear (lrop);
   }
}

/*****************************************************************************/

void mpfpx_div_qr (mpfpx_t q, mpfpx_t r, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv)

{
   mpfpx_t lq, lr;
   mpfp_t  inv, tmp_mpfp;
      /* the inverse of the leading coefficient */
   int     i, j;
   int    monic = false;
      /* true if op2 is monic */

   mpfpx_init (lq);
   mpfpx_init_set (lr, op1);
   mpfp_init (inv);
   mpfp_init (tmp_mpfp);

   lq->deg = lr->deg - op2->deg;
   if (lq->deg < 0)
      mpfpx_deg_set (lq, -1);
   else if (lc_inv != NULL)
      mpfp_set (inv, lc_inv);
   else
      mpfp_inv (inv, op2->coeff [op2->deg]);
   if (mpfp_cmp_ui (inv, 1ul) == 0)
      monic = true;

   if (__mpfpx_arithmetic == MPFPX_KARATSUBA)
   /* It would not hurt to use this part of code even with naive arithmetic. */
   /* Then mpfpx_mul_low also uses naive arithmetic, and the number of       */
   /* operations remains unchanged altogether.                               */
   {
      /* first compute the quotient */
      mpfpx_div_q (lq, op1, op2, inv);
      /* then compute the remainder as lr = op1 - lq*op2 */
      mpfpx_mul_low (lr, lq, op2, op2->deg - 1);
      mpfpx_sub (lr, op1, lr);
      lr->deg = op2->deg - 1;
   }
   else
   for (i = lq->deg; i >= 0; i--)
   {
      /* compute the coefficient i of q */
      if (monic)
         mpfp_set (lq->coeff [i], lr->coeff [lr->deg]);
      else
         mpfp_mul (lq->coeff [i], lr->coeff [lr->deg], inv);
      /* multiply back */
      lr->deg --;
      for (j = 0; j < op2->deg; j++)
      {
         mpfp_mul (tmp_mpfp, lq->coeff [i], op2->coeff [j]);
         mpfp_sub (lr->coeff [i+j], lr->coeff [i+j], tmp_mpfp);
      }
   }
   mpfpx_normalise (lr);

   mpfpx_set (q, lq);
   mpfpx_set (r, lr);

   mpfpx_clear (lq);
   mpfpx_clear (lr);
   mpfp_clear (inv);
   mpfp_clear (tmp_mpfp);
}

/*****************************************************************************/

void mpfpx_div_qr_newton (mpfpx_t q, mpfpx_t r, mpfpx_t op1, mpfpx_t op2)

{
   mpfpx_t tmp, lop1, lop2;
   int     q_deg = op1->deg - op2->deg;

   mpfpx_init (tmp);
   mpfpx_init_set (lop1, op1);
   mpfpx_init_set (lop2, op2);

   if (op1->deg < op2->deg)
   {
      q->deg = -1;
      mpfpx_set (r, op1);
   }
   else
   {
      printf ("a "); mpfpx_out (op1);
      printf ("b "); mpfpx_out (op2);
      mpfpx_rev (tmp, op2, op2->deg);
      printf ("rev b "); mpfpx_out (tmp);
      mpfpx_newton (tmp, tmp, q_deg + 1);
      printf ("inv "); mpfpx_out (tmp);
      mpfpx_rev (q, op1, op1->deg);
      q->deg = MIN(q_deg, q->deg);
      mpfpx_mul (q, q, tmp);
      q->deg = MIN(q_deg, q->deg);
      mpfpx_rev (q, q, q_deg);
      lop2->deg--;
      mpfpx_mul (tmp, q, lop2);
      mpfpx_sub (r, lop1, tmp);
      r->deg = lop2->deg;
      lop2->deg++;
      mpfpx_normalise (r);
   }

   mpfpx_clear (tmp);
   mpfpx_clear (lop1);
   mpfpx_clear (lop2);
}

/*****************************************************************************/

void mpfpx_div_q (mpfpx_t q, mpfpx_t op1, mpfpx_t op2, mpfp_t lc_inv)

{
   mpfpx_t lq, lr;
   mpfp_t  inv, tmp_mpfp;
      /* the inverse of the leading coefficient */
   int     i, j;
   int    monic = false;
      /* true if op2 is monic */

   mpfpx_init (lq);
   mpfpx_init_set (lr, op1);
   mpfp_init (inv);
   mpfp_init (tmp_mpfp);

   mpfpx_deg_set (lq, lr->deg - op2->deg);
   if (lq->deg < 0)
      mpfpx_deg_set (lq, -1);
   else if (lc_inv != NULL)
      mpfp_set (inv, lc_inv);
   else
      mpfp_inv (inv, op2->coeff [op2->deg]);
   if (mpfp_cmp_ui (inv, 1ul) == 0)
      monic = true;

   for (i = lq->deg; i >= 0; i--)
   {
      /* compute the coefficient i of q */
      if (monic)
         mpfp_set (lq->coeff [i], lr->coeff [lr->deg]);
      else
         mpfp_mul (lq->coeff [i], lr->coeff [lr->deg], inv);
      /* multiply back */
      lr->deg --;
      for (j = MAX(op2->deg - i, 0); j < op2->deg; j++)
      {
         mpfp_mul (tmp_mpfp, lq->coeff [i], op2->coeff [j]);
         mpfp_sub (lr->coeff [i+j], lr->coeff [i+j], tmp_mpfp);
      }
   }

   mpfpx_set (q, lq);

   mpfpx_clear (lq);
   mpfpx_clear (lr);
   mpfp_clear (inv);
   mpfp_clear (tmp_mpfp);
}

/*****************************************************************************/

void mpfpx_div_r (mpfpx_t r, mpfpx_t op1, mpfpx_t op2, mpfp_t lc_inv)
   /* This is a verbatim copy of mpfpx_div_qr, except for the assignment of q. */

{
   mpfpx_t lq, lr;
   mpfp_t  inv, tmp_mpfp;
      /* the inverse of the leading coefficient */
   int     i, j;
   int    monic = false;
      /* true if op2 is monic */

   mpfpx_init (lq);
   mpfpx_init_set (lr, op1);
   mpfp_init (inv);
   mpfp_init (tmp_mpfp);

   lq->deg = lr->deg - op2->deg;
   if (lq->deg < 0)
      mpfpx_deg_set (lq, -1);
   else if (lc_inv != NULL)
      mpfp_set (inv, lc_inv);
   else
      mpfp_inv (inv, op2->coeff [op2->deg]);
   if (mpfp_cmp_ui (inv, 1ul) == 0)
      monic = true;

   if (__mpfpx_arithmetic == MPFPX_KARATSUBA)
   /* It would not hurt to use this part of code even with naive arithmetic. */
   /* Then mpfpx_mul_low also uses naive arithmetic, and the number of       */
   /* operations remains unchanged altogether.                               */
   {
      /* first compute the quotient */
      mpfpx_div_q (lq, op1, op2, inv);
      /* then compute the remainder as lr = op1 - lq*op2 */
      mpfpx_mul_low (lr, lq, op2, op2->deg - 1);
      mpfpx_sub (lr, op1, lr);
      lr->deg = op2->deg - 1;
   }
   else
   for (i = lq->deg; i >= 0; i--)
   {
      /* compute the coefficient i of q */
      if (monic)
         mpfp_set (lq->coeff [i], lr->coeff [lr->deg]);
      else
         mpfp_mul (lq->coeff [i], lr->coeff [lr->deg], inv);
      /* multiply back */
      lr->deg --;
      for (j = 0; j < op2->deg; j++)
      {
         mpfp_mul (tmp_mpfp, lq->coeff [i], op2->coeff [j]);
         mpfp_sub (lr->coeff [i+j], lr->coeff [i+j], tmp_mpfp);
      }
   }
   mpfpx_normalise (lr);

   mpfpx_set (r, lr);

   mpfpx_clear (lq);
   mpfpx_clear (lr);
   mpfp_clear (inv);
   mpfp_clear (tmp_mpfp);
}

/*****************************************************************************/

void mpfpx_divexact_jebelean (mpfpx_t q, mpfpx_t op1, mpfpx_t op2)

{
   mpfpx_t lop1, lop2;
   mpfp_t  inv, tmp_mpfp;
      /* the inverse of the constant coefficient */
   int     i, j;

   mpfpx_init_set (lop1, op1);
   mpfpx_init_set (lop2, op2);
   mpfp_init (inv);
   mpfp_init (tmp_mpfp);

   mpfpx_deg_set (q, lop1->deg - lop2->deg);
   if (q->deg < 0)
      mpfpx_deg_set (q, -1);
   else
      mpfp_inv (inv, lop2->coeff [0]);

   for (i = 0; i <= q->deg; i++)
   {
      /* compute the coefficient i of q */
      mpfp_mul (q->coeff [i], lop1->coeff [i], inv);
      /* multiply back */
      for (j = 1; j <= q->deg - i; j++)
      {
         mpfp_mul (tmp_mpfp, q->coeff [i], lop2->coeff [j]);
         mpfp_sub (lop1->coeff [i+j], lop1->coeff [i+j], tmp_mpfp);
      }
   }

   mpfpx_clear (lop1);
   mpfpx_clear (lop2);
   mpfp_clear (inv);
   mpfp_clear (tmp_mpfp);
}

/*****************************************************************************/

void mpfpx_divexact (mpfpx_t q, mpfpx_t op1, mpfpx_t op2, mpfp_t lc_inv)

{
   mpfpx_div_q (q, op1, op2, lc_inv);
}

/*****************************************************************************/

void mpfpx_neg (mpfpx_t rop, mpfpx_t op)

{
   int i;

   rop->deg = op->deg;
   for (i = 0; i <= op->deg; i++)
      mpfp_neg (rop->coeff [i], op->coeff [i]);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfpx_gcd_ext (mpfpx_t d, mpfpx_t u, mpfpx_t v, mpfpx_t a, mpfpx_t b,
   int steps, mpfp_t lc_inv)

{
   mpfpx_t la, q, u0, v0, tmp;
   int    found = false;

   mpfpx_init_set (la, a);
   mpfpx_set (d, b);
   mpfpx_init (q);
   if (u != NULL)
   {
      mpfpx_init_set_ui (u0, 1ul);
      mpfpx_set_ui (u, 0ul);
   }
   if (v != NULL)
   {
      mpfpx_init_set_ui (v0, 0ul);
      mpfpx_set_ui (v, 1ul);
   }
   /* we always have la = u0 * a + v0 * b and d = u * a + v * b */
   mpfpx_init (tmp);

   while (!found && d->deg > 0 && steps > 0)
   {
      /* compute remainder */
      if (d->deg == b->deg)
         mpfpx_div_qr (q, tmp, la, d, lc_inv);
      else
         mpfpx_div_qr (q, tmp, la, d, NULL);
      steps--;
      if (tmp->deg >= 0)
      {
         mpfpx_set (la, d);
         mpfpx_set (d, tmp);
         if (u != NULL)
         {
            /* compute u */
            mpfpx_set (tmp, u0);
            mpfpx_set (u0, u);
            mpfpx_set (u, tmp);
            mpfpx_mul (tmp, q, u0);
            mpfpx_sub (u, u, tmp);
         }
         if (v != NULL)
         {
            /* compute v */
            mpfpx_set (tmp, v0);
            mpfpx_set (v0, v);
            mpfpx_set (v, tmp);
            mpfpx_mul (tmp, q, v0);
            mpfpx_sub (v, v, tmp);
         }
      }
      else
         found = true;
   }

   mpfpx_clear (la);
   mpfpx_clear (q);
   if (u != NULL)
   {
      mpfpx_clear (u0);
      mpfpx_normalise (u);
   }
   if (v != NULL)
   {
      mpfpx_clear (v0);
      mpfpx_normalise (v);
   }
   mpfpx_clear (tmp);
}

/*****************************************************************************/

/*
void mpfpx_gcd_ext (mpfpx_t d, mpfpx_t u, mpfpx_t v, mpfpx_t a, mpfpx_t b,
   int steps, mpfp_t lc_inv)
   // normalising all intermediate remainders as in GG, Algorithm 3.6;
   // correct, but very inefficient

{
   mpfp_t  one, inv;
   mpfpx_t la, q, u0, v0, tmp;
   bool    found = false;

   mpfp_init_set_ui (one, 1);
   mpfp_init (inv);
   mpfpx_init (la);
   mpfpx_init (q);

   mpfp_inv (inv, a->coeff [a->deg]);
   mpfpx_monicise (la, a, inv);
   if (u != NULL)
   {
      mpfpx_init_set_fp (u0, inv);
      mpfpx_set_ui (u, 0);
   }
   if (lc_inv == NULL)
      mpfp_inv (inv, b->coeff [b->deg]);
   else
      mpfp_set (inv, lc_inv);
   mpfpx_monicise (d, b, inv);
   if (v != NULL)
   {
      mpfpx_init_set_ui (v0, 0);
      mpfpx_set_fp (v, inv);
   }
   // we always have la = u0 * a + v0 * b and d = u * a + v * b, la and d monic
   mpfpx_init (tmp);

   while (!found && d->deg > 0 && steps > 0)
   {
      // compute remainder
      mpfpx_div_qr (q, tmp, la, d, one);
      steps--;
      if (tmp->deg >= 0)
      {
         mpfp_inv (inv, tmp->coeff [tmp->deg]);
         mpfpx_set (la, d);
         mpfpx_monicise (d, tmp, inv);
         if (u != NULL)
         {
            // compute u
            mpfpx_set (tmp, u0);
            mpfpx_set (u0, u);
            mpfpx_set (u, tmp);
            mpfpx_mul (tmp, q, u0);
            mpfpx_sub (u, u, tmp);
            mpfpx_set_fp (tmp, inv);
            mpfpx_mul (u, u, tmp);
         }
         if (v != NULL)
         {
            // compute v
            mpfpx_set (tmp, v0);
            mpfpx_set (v0, v);
            mpfpx_set (v, tmp);
            mpfpx_mul (tmp, q, v0);
            mpfpx_sub (v, v, tmp);
            mpfpx_set_fp (tmp, inv);
            mpfpx_mul (v, v, tmp);
         }
      }
      else
         found = true;
   }

   mpfp_clear (one);
   mpfp_clear (inv);
   mpfpx_clear (la);
   mpfpx_clear (q);
   if (u != NULL)
      mpfpx_clear (u0);
   if (v != NULL)
      mpfpx_clear (v0);
   mpfpx_clear (tmp);
}
*/

/*****************************************************************************/

void mpfpx_invert (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, mpfp_t lc_inv)

{
   mpfpx_t d;
   int     i;

   mpfpx_init (d);

   mpfpx_gcd_ext (d, rop, NULL, op1, op2, MPFPX_INFTY, lc_inv);
   if (d->deg != 0)
   {
      printf ("*** Error in mpfpx_invert!\n");
      exit (1);
   }

   mpfp_inv (d->coeff [0], d->coeff [0]);
   for (i = 0; i <= rop->deg; i++)
      mpfp_mul (rop->coeff [i], rop->coeff [i], d->coeff [0]);

   mpfpx_clear (d);
}

/*****************************************************************************/
