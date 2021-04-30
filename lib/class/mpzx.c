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

