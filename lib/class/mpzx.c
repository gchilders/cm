/*

mpzx.c - code for handling polynomials with mpz coefficients
Unless stated otherwise, the functions and their prototypes are
inspired by mpfrcx.

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

#include "cm_class-impl.h"

/*****************************************************************************/

static void quadratic_basis (ctype omega, int_cl_t d);

/*****************************************************************************/

static void quadratic_basis (ctype omega, int_cl_t d)
   /* Compute the standard omega such that [1, omega] is an integral basis
      of the imaginary quadratic order with discriminant d, that is,
      omega = sqrt(d)/2 if d is divisible by 4,
      omega = (1+sqrt(d))/2 otherwise. */
{
   fsqrt_ui (cimagref (omega), (unsigned long int) (-d));
   if (cm_classgroup_mod (d, (uint_cl_t) 4) == 0)
      fset_ui (crealref (omega), 0);
   else
      fset_ui (crealref (omega), 1);
   cdiv_2ui (omega, omega, 1);
}

/*****************************************************************************/

bool cm_nt_cget_quadratic (mpz_t out1, mpz_t out2, ctype in, int_cl_t d)
{
   ctype omega;
   bool ok;

   cinit (omega, cget_prec (in));
   quadratic_basis (omega, d);

   ok = cm_nt_cget_zz (out1, out2, in, omega);

   return ok;
}

/*****************************************************************************/

void mpzx_init (mpzx_ptr f, int deg)
   /* Initialise a polynomial of degree deg, that is, with deg+1
      coefficients, and also set its degree. Unlike for mpfrx, for
      instance, there is no distinction between the available space
      and the degree, since in our applications the degree is usually
      known beforehand. */
{
   int i;

   f->deg = deg;
   f->coeff = (mpz_t *) malloc ((deg + 1) * sizeof (mpz_t));
   for (i = 0; i <= deg; i++)
      mpz_init (f->coeff [i]);
}

/*****************************************************************************/

void mpzx_clear (mpzx_ptr f)
{
   int i;

   for (i = 0; i <= f->deg; i++)
      mpz_clear (f->coeff [i]);
   free (f->coeff);
}

/*****************************************************************************/

size_t mpzx_out_str (FILE* stream, int base, mpzx_srcptr f)
{
   if (stream == NULL)
      stream = stdout;

   if (f->deg == -1) {
      fprintf (stream, "(-1 0)");
      return (size_t) 6u;
   }
   else {
      size_t s = f->deg + 3u; /* for '(', ')' and the spaces */
      int i;

      fprintf (stream, "(");
      s += fprintf (stream, "%i ", f->deg);
      for (i = f->deg; i >= 0; i--) {
         s += mpz_out_str (stream, base, f->coeff [i]);
         if (i > 0)
            fprintf (stream, " ");
         else
            fprintf (stream, ")");
      }

      return s;
   }
}

/*****************************************************************************/

void mpzx_print_pari (FILE* file, mpzx_srcptr f, char *var)
   /* Print the polynomial f in a format understood by PARI.
      var contains the variable name; if it is NULL, the function uses "x". */
{
   const char *x = (var == NULL ? "x" : var);
   int i, cmp;

   if (f->deg == -1)
      fprintf (file, "0");
   else {
      for (i = f->deg; i >= 0; i--) {
         cmp = mpz_cmp_ui (f->coeff [i], 0);
         if (cmp != 0) {
            if (cmp > 0 && i != f->deg)
               fprintf (file, "+");
            if (mpz_cmp_ui (f->coeff [i], 1) == 0) {
               if (i == 0)
                  fprintf (file, "1");
            }
            else if (mpz_cmp_si (f->coeff [i], -1) == 0)
               if (i == 0)
                  fprintf (file, "-1");
               else
                  fprintf (file, "-");
            else
               mpz_out_str (file, 10, f->coeff [i]);
            if (i > 0 && mpz_cmpabs_ui (f->coeff [i], 1) != 0)
               fprintf (file, "*");
            if (i >= 1)
               fprintf (file, "%s", x);
            if (i >= 2)
               fprintf (file, "^%i", i);
         }
      }
   }
}


/*****************************************************************************/

void mpzx_tower_init (mpzx_tower_ptr twr, int levels, int *d)
{
   int i, j, deg;

   twr->levels = levels;
   twr->d = (int *) malloc (levels * sizeof (int));
   deg = 1;
   for (i = 0; i < levels; i++) {
      twr->d [i] = d [i];
      deg *= d [i];
   }
   twr->deg = deg;
   twr->W = (mpzx_t **) malloc (levels * sizeof (mpzx_t *));
   twr->W [0] = (mpzx_t *) malloc (1 * sizeof (mpzx_t));
   mpzx_init (twr->W [0][0], d [0]);
   deg = 1;
   for (i = 1; i < levels; i++) {
      twr->W [i] = (mpzx_t *) malloc ((d [i] + 1) * sizeof (mpzx_t));
      deg *= d [i - 1];
      for (j = 0; j <= d [i]; j++)
         mpzx_init (twr->W [i][j], deg-1);
   }
}

/*****************************************************************************/

void mpzx_tower_clear (mpzx_tower_ptr twr)
{
   int i, j;

   mpzx_clear (twr->W [0][0]);
   free (twr->W [0]);
   for (i = 1; i < twr->levels; i++) {
      for (j = 0; j <= twr->d [i]; j++)
         mpzx_clear (twr->W [i][j]);
      free (twr->W [i]);
   }
   free (twr->W);
   free (twr->d);
}

/*****************************************************************************/

void mpzx_tower_print_pari (FILE* file, mpzx_tower_srcptr twr, char *fun,
   char *var)
   /* Print the number field tower twr in a format understood by PARI.
      fun contains the base name used for the polynomials; if it is NULL,
      the function uses "f". var contains the base name for the variables;
      if it is NULL, the function uses "x". */

{
   const char *f = (fun == NULL ? "f" : fun);
   const char *x = (var == NULL ? "x" : var);
   int i, j;
   char xi [22]; /* long enough to hold x18446744073709551615; assumes that
                    var contains only one character */

   fprintf (file, "%s1 = ", f);
   sprintf (xi, "%s1", x);
   mpzx_print_pari (file, twr->W [0][0], xi);
   printf (";\n");
   for (i = 1; i < twr->levels; i++) {
      printf ("%s%i = ", f, i+1);
      sprintf (xi, "%s%u", x, (unsigned int) i);
      for (j = twr->d [i]; j >= 0; j--) {
         if (j < twr->d [i])
            printf ("+");
         printf ("(");
         mpzx_print_pari (file, twr->W [i][j], xi);
         if (j >= 2)
            printf (")*%s%i^%i", x, i+1, j);
         else if (j == 1)
            printf (")*%s%i", x, i+1);
         else
            printf (")");
      }
      printf (";\n");
   }
}

/*****************************************************************************/

