/*

modpol.c - code for handling modular polynomials

Copyright (C) 2009, 2012, 2021 Andreas Enge

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

static unsigned long int read_gz_ui (gzFile f);
static void read_gz_mpz (mpz_t rop, gzFile f);

/*****************************************************************************/

static unsigned long int read_gz_ui (gzFile f)
   /* reads and returns an unsigned long integer in base 10 from the gzipped */
   /* file f                                                                 */
{
   char c;
   unsigned long int res = 0ul;

   do {
      c = gzgetc (f);
      if (c == EOF) {
         printf ("*** EOF in 'modpol_read_gz_ui'\n");
         exit (1);
      }
   } while (isspace (c));
   do {
      if (c < '0' || c > '9') {
         printf ("*** Non-digit in 'modpol_read_gz_ui'\n");
         exit (1);
      }
      res = 10 * res + (c - '0');
      c = gzgetc (f);
   } while (!(isspace (c) || c == EOF));

   return res;
}

/*****************************************************************************/

static void read_gz_mpz (mpz_t rop, gzFile f)
   /* reads a coefficient of type mpz in base 10 from the gzipped file f     */
{
   char c;
   int sign = 1;

   mpz_set_ui (rop, 0);

   do {
      c = gzgetc (f);
      if (c == EOF) {
         printf ("*** EOF in 'modpol_read_gz_mpz'\n");
         exit (1);
      }
   } while (isspace (c));
   if (c == '-') {
      sign = -1;
      c = gzgetc (f);
      if (c == EOF) {
         printf ("*** EOF in 'modpol_read_gz_mpz'\n");
         exit (1);
      }
   }
   do {
      if (c < '0' || c > '9') {
         printf ("*** Non-digit in 'modpol_read_gz_mpz'\n");
         exit (1);
      }
      mpz_mul_ui (rop, rop, 10ul);
      mpz_add_ui (rop, rop, (unsigned long int) (c - '0'));
      c = gzgetc (f);
   } while (!(isspace (c) || c == EOF));

   if (sign == -1)
      mpz_neg (rop, rop);
}

/*****************************************************************************/

mpz_t* cm_modpol_read_specialised_mod (int* n, int level, char type, mpz_t p,
   mpz_t x, const char* datadir)
   /* returns the modular polynomial of the given type and level modulo p    */
   /* and evaluated in the first argument x; n is replaced by its degree  */

{
   char filename [255];
   gzFile f;
   int lev, N, i, k;
   char c;
   mpz_t* xpow;
   mpz_t tmp;
   mpz_t* res;

   sprintf (filename, "%s/%cf/%cf_%.4i.dat.gz", datadir,
      type, type, level);
   cm_file_gzopen_read (&f, filename);

   lev = read_gz_ui (f);
   if (lev != level) {
      printf ("*** Trying to read modular polynomial of level %i ", level);
      printf ("in a file for the level %i!\n", lev);
      exit (1);
   }
   c = gzgetc (f);
   if (c != type) {
      printf ("*** Trying to read modular polynomial of type %c ", type);
      printf ("in a file for the type %c!\n", c);
      exit (1);
   }

   N = read_gz_ui (f);
   xpow = (mpz_t *) malloc ((N + 1) * sizeof (mpz_t));
   mpz_init_set_ui (xpow [0], 1ul);
   for (i = 1; i <= N; i++) {
      mpz_init (xpow [i]);
      mpz_mul (xpow [i], xpow [i-1], x);
      mpz_mod (xpow [i], xpow [i], p);
   }
   mpz_init (tmp);

   *n = read_gz_ui (f);
   res = (mpz_t*) malloc ((*n+1) * sizeof (mpz_t));
   for (i = 0; i <= *n; i++)
      mpz_init_set_ui (res [i], 0ul);
   do {
      i = read_gz_ui (f);
      k = read_gz_ui (f);
      read_gz_mpz (tmp, f);
      mpz_mul (tmp, tmp, xpow [i]);
      mpz_add (res [k], res [k], tmp);
      mpz_mod (res [k], res [k], p);
   } while (k != 0 || i != 0);
      /* we assume that the last entry in the file is the constant one */

   for (i = 0; i <= N; i++)
      mpz_clear (xpow [i]);
   free (xpow);
   mpz_clear (tmp);

   cm_file_gzclose (f);

   return res;
}

/*****************************************************************************/

void cm_modpol_print_pari (int level, char type, const char* datadir)
      /* prints the modular polynomial of the given type and level in the    */
      /* pari seadata format                                                 */

{
   char filename [255];
   gzFile f;
   int lev, i_old, i, k;
   char c;
   mpz_t tmp;

   if (type != 'a') {
      printf ("*** Trying to read modular polynomial of type %c ", type);
      printf ("instead of 'a'!\n");
      exit (1);
   }
   sprintf (filename, "%s/%cf/%cf_%.4i.dat.gz", datadir,
            type, type, level);
   cm_file_gzopen_read (&f, filename);

   lev = read_gz_ui (f);
   if (lev != level) {
      printf ("*** Trying to read modular polynomial of level %i ", level);
      printf ("in a file for the level %i!\n", lev);
      exit (1);
   }
   c = gzgetc (f);
   if (c != type) {
      printf ("*** Trying to read modular polynomial of type '%c' ", type);
      printf ("in a file for the type %c!\n", c);
      exit (1);
   }

   /* skip N and n */
   read_gz_ui (f);
   read_gz_ui (f);
   mpz_init (tmp);

   printf ("[%i, \"A\", [", level);
   i_old = level + 2;
   do {
      i = read_gz_ui (f);
      k = read_gz_ui (f);
      read_gz_mpz (tmp, f);
      if (i != i_old && k != 0)
         printf ("[");
      mpz_out_str (stdout, 10, tmp);
      if (k != 0)
         printf (", ");
      else {
         if (i == i_old)
            printf ("]");
         if (i != 0)
            printf (", ");
      }
      i_old = i;
   } while (k != 0 || i != 0);
   /* we assume that the last entry in the file is the constant one */
   printf ("]]");

   cm_file_gzclose (f);
}

/*****************************************************************************/

void cm_modpol_print_magma (int level, char type, const char* datadir)
      /* prints the modular polynomial of the given type and level in magma  */
      /* format                                                              */

{
   char filename [255];
   gzFile f;
   int lev, i, k;
   char c;
   mpz_t tmp;

   sprintf (filename, "%s/%cf/%cf_%.4i.dat.gz", datadir,
            type, type, level);
   cm_file_gzopen_read (&f, filename);

   lev = read_gz_ui (f);
   if (lev != level) {
      printf ("*** Trying to read modular polynomial of level %i ", level);
      printf ("in a file for the level %i!\n", lev);
      exit (1);
   }
   c = gzgetc (f);
   if (c != type) {
      printf ("*** Trying to read modular polynomial of type %c ", type);
      printf ("in a file for the type %c!\n", c);
      exit (1);
   }

   /* skip N and n */
   read_gz_ui (f);
   read_gz_ui (f);
   mpz_init (tmp);

   printf ("phi :=");
   do {
      i = read_gz_ui (f);
      k = read_gz_ui (f);
      read_gz_mpz (tmp, f);
      printf (" +(");
      mpz_out_str (stdout, 10, tmp);
      printf (")*F^%i*J^%i", i, k);
   } while (k != 0 || i != 0);
   /* we assume that the last entry in the file is the constant one */
   printf (";\n");

   cm_file_gzclose (f);
}

/*****************************************************************************/
