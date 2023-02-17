/*

mpzx.c - code for handling polynomials with mpz coefficients
Unless stated otherwise, the functions and their prototypes are
inspired by mpfrcx.

Copyright (C) 2021, 2022, 2023 Andreas Enge

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

/*****************************************************************************/

static void quadratic_basis (ctype omega, int_cl_t d);
static bool mpcx_get_mpzxx (mpzx_ptr g, mpzx_ptr h, mpcx_srcptr f,
   ctype omega);
static uint64_t mpz_hash (mpz_srcptr z);
static uint64_t mpzx_hash (mpzx_srcptr f);

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

void mpzx_set_deg (mpzx_ptr f, int deg)
   /* Modify the degree of f to deg; this destroys the content of f if
      the old and the new degree are different. Otherwise the function
      does nothing. */
{
   if (deg != f->deg) {
      mpzx_clear (f);
      mpzx_init (f, deg);
   }
}

/*****************************************************************************/

void mpzx_set (mpzx_ptr f, mpzx_srcptr g)
   /* Set f to g. */
{
   int i;

   mpzx_set_deg (f, g->deg);
   for (i = 0; i <= g->deg; i++)
      mpz_set (f->coeff [i], g->coeff [i]);
}

/*****************************************************************************/

int mpzx_cmp (mpzx_srcptr f, mpzx_srcptr g)
   /* Compare f and g. The return value is 0 when both are equal, -1 when
      f is smaller than g and +1 when f is larger than g for some arbitrary
      ordering. */
{
   int i, cmp;

   if (f->deg < g->deg)
      return -1;
   else if (f->deg > g->deg)
      return 1;
   else /* the degrees are equal */ {
      for (i = f->deg; i >= 0; i--) {
         cmp = mpz_cmp (f->coeff [i], g->coeff [i]);
         if (cmp != 0)
            return cmp;
      }
      return 0;
   }
}

/*****************************************************************************/

void mpzx_mod (mpzx_ptr g, mpzx_srcptr f, mpz_srcptr p)
   /* Compute g = f modulo p. */
{
   int i;

   mpzx_set_deg (g, f->deg);
   for (i = 0; i <= f->deg; i++)
      mpz_mod (g->coeff [i], f->coeff [i], p);
}

/*****************************************************************************/

bool cm_mpfrx_get_mpzx (mpzx_ptr g, mpfrx_srcptr f)
   /* Try to round f to g; the return value reflects the success of the
      operation. It is assumed that g already has the same degree as f. */
{
   int k;
   bool ok = true;

   for (k = 0; k <= mpfrx_get_deg (f); k++)
      ok &= cm_nt_fget_z (g->coeff [k], mpfrx_get_coeff (f, k));

   return ok;
}

/*****************************************************************************/

static bool mpcx_get_mpzxx (mpzx_ptr g, mpzx_ptr h, mpcx_srcptr f,
   ctype omega)
   /* Try to round f such that f = g + omega*h. It is assumed that g and h
      already have the same degree as f. The return value reflects the
      success of the operation. */
{
   int k;
   bool ok = true;

   for (k = 0; k <= mpcx_get_deg (f); k++)
      ok &= cm_nt_cget_zz (g->coeff [k], h->coeff [k],
         mpcx_get_coeff (f, k), omega);

   return ok;
}

/*****************************************************************************/

bool cm_mpcx_get_quadraticx (mpzx_ptr g, mpzx_ptr h, mpcx_srcptr f,
   int_cl_t d)
   /* Try to round f such that f = g + omega*h, where omega is the second
      element of the standard basis attached to the fundamental
      discriminant d. It is assumed that g and h already have the same
      degree as f. The return value reflects the success of the operation. */
{
   ctype omega;
   bool ok;

   cinit (omega, cget_prec (mpcx_get_coeff (f, 0)));
   quadratic_basis (omega, d);

   ok = mpcx_get_mpzxx (g, h, f, omega);

   cclear (omega);

   return ok;
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
      var contains the variable name; if it is NULL, the function
      uses "x". */
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

void mpzxx_print_pari (FILE* file, mpzx_srcptr g, mpzx_srcptr h, char *var)
   /* Print the polynomial f = g + omega*h in a format understood by PARI.
      var contains the variable name; if it is NULL, the function
      uses "x". */
{
   fprintf (file, "(");
   mpzx_print_pari (file, g, var);
   fprintf (file, ")+o*(");
   mpzx_print_pari (file, h, var);
   fprintf (file, ")");
}

/*****************************************************************************/
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

bool cm_mpfrx_tower_get_mpzx_tower (mpzx_tower_ptr tz, mpfrx_tower_srcptr tf)
   /* Try to round the floating point field tower tf to the integral tower
      tz, assuming that tz is already initialised with the correct degrees.
      The return value reflects the success of the operation. */
{
   int i, j;
   bool ok;

   ok = cm_mpfrx_get_mpzx (tz->W [0][0], tf->W [0][0]);
   for (i = 1; i < tf->levels; i++)
      for (j = tf->d [i]; j >= 0; j--)
         ok &= cm_mpfrx_get_mpzx (tz->W [i][j], tf->W [i][j]);

   return ok;
}

/*****************************************************************************/
bool cm_mpcx_tower_get_quadratic_tower (mpzx_tower_ptr t1,
   mpzx_tower_ptr t2, mpcx_tower_srcptr tc, int_cl_t d)
   /* Try to round the complex floating point field tower tc as
      "tc = t1 + omega*t2", where omega is the second basis element of the
      standard basis for the imaginary-quadratic number field of
      discriminant d. t1 and t2 are assumed to be initialised with the
      correct degree sequence. The return value reflects the success of
      the operation. */
{
   int i, j;
   ctype omega;
   bool ok;

   cinit (omega, cget_prec (mpcx_get_coeff (tc->W [0][0], 0)));
   quadratic_basis (omega, d);

   ok = cm_mpcx_get_quadraticx (t1->W [0][0], t2->W [0][0], tc->W [0][0], d);
   for (i = 1; i < tc->levels; i++)
      for (j = tc->d [i]; j >= 0; j--)
         ok = cm_mpcx_get_quadraticx (t1->W [i][j], t2->W [i][j], tc->W [i][j], d);

   cclear (omega);

   return ok;

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

void mpzxx_tower_print_pari (FILE* file, mpzx_tower_srcptr g,
   mpzx_tower_srcptr h, char *fun,
   char *var)
   /* Print the number field tower f = g + omega * h in a format
      understood by PARI.
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
   mpzxx_print_pari (file, g->W [0][0], h->W [0][0], xi);
   printf (";\n");
   for (i = 1; i < g->levels; i++) {
      printf ("%s%i = ", f, i+1);
      sprintf (xi, "%s%u", x, (unsigned int) i);
      for (j = g->d [i]; j >= 0; j--) {
         if (j < g->d [i])
            printf ("+");
         printf ("(");
         mpzxx_print_pari (file, g->W [i][j], h->W [i][j], xi);
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
/*****************************************************************************/

static uint64_t mpz_hash (mpz_srcptr z)
   /* Return a hash value for z, without special thought to potential
      collisions. */
{
   int len;
   mp_limb_t h;
   int i;

   len = (z->_mp_size >= 0 ? z->_mp_size : -z->_mp_size);
   h = (mp_limb_t) z->_mp_size;
   for (i = 0; i < len; i++)
      h ^= z->_mp_d [i];

   return (uint64_t) h;
}

/*****************************************************************************/

static uint64_t mpzx_hash (mpzx_srcptr f)
   /* Return a hash value for f, without special thought to potential
      collisions. */
{
   uint64_t h, h0;
   unsigned int deg, shift, i;

   h = (uint64_t) f->deg;
   if (f->deg >= 0) {
      deg = (unsigned int) f->deg;
      for (i = 0; i <= deg; i++) {
         h0 = mpz_hash (f->coeff [i]);
         /* Circular shift by the position to break symmetries, so that
            at least polynomials such as X^2+X-1 and X^2-X+1 do not get
            the same hash value. */
         shift = i % 64;
         if (shift != 0)
            h0 = (h0 << shift) | (h0 >> (64 - shift));
         h ^= h0;
      }
   }

   return h;
}

/*****************************************************************************/

uint64_t mpzx_mod_hash (mpzx_srcptr f, mpz_srcptr p)
   /* Returns a hash value for f modulo p, without special thought to
      potential collisions. */
{
   return mpzx_hash (f) ^ mpz_hash (p);
}

/*****************************************************************************/

