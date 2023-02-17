/*

file.c - code for handling (gzipped) files

Copyright (C) 2009, 2012, 2021, 2022, 2023 Andreas Enge

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

static bool read_ecpp_cert1_line (FILE *f, mpz_t *line, cm_stat_t stat);
static int read_ecpp_cert2_line (FILE *f, mpz_t *line, cm_stat_t stat);
static bool write_stat (FILE *f, cm_stat_t stat);
static bool read_stat (FILE *f, cm_stat_t stat);
static void mpz_out_hex (FILE *f, mpz_t z);

/*****************************************************************************/
/*                                                                           */
/* Printing data augmented with the MPI rank and wallclock time.             */
/*                                                                           */
/*****************************************************************************/

void cm_file_printf (char *fmt, ...)
{
   va_list ap;

#ifdef WITH_MPI
   int rank;
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#endif

   va_start (ap, fmt);
#ifndef WITH_MPI
   vprintf (fmt, ap);
#else
   printf ("[%4i] (%li)  ", rank, time (NULL) - cm_mpi_zero);
   vprintf (fmt, ap);
#endif
   va_end (ap);
   fflush (stdout);
}

/*****************************************************************************/

bool cm_file_open_write (FILE **f, const char *filename)
{
   bool res;

   *f = fopen (filename, "w");
   if (*f == NULL) {
      printf ("Could not open file '%s' for writing.\n", filename);
      res = false;
   }
   else {
      printf ("Writing to '%s'.\n", filename);
      res = true;
   }
   fflush (stdout);

   return res;
}

/*****************************************************************************/

bool cm_file_open_read (FILE **f, const char *filename)
{
   bool res;

   *f = fopen (filename, "r");
   if (*f == NULL) {
      printf ("Could not open file '%s' for reading.\n", filename);
      res = false;
   }
   else {
      printf ("Reading from '%s'.\n", filename);
      res = true;
   }
   fflush (stdout);

   return res;
}

/*****************************************************************************/

bool cm_file_open_read_write (FILE **f, const char *filename)
{
   bool res;
   *f = fopen (filename, "r+");
   if (*f == NULL) {
      printf ("Could not open file '%s' for reading.\n", filename);
      res = false;
   }
   else
      res = true;
   fflush (stdout);

   return res;
}

/*****************************************************************************/

void cm_file_close (FILE *f)
{
   fclose (f);
}

/*****************************************************************************/

void cm_file_gzopen_write (gzFile *f, const char *filename)
{
   *f = gzopen (filename, "w9");
   if (*f == NULL) {
      printf ("Could not open file '%s' for writing.\n", filename);
      exit (1);
   }
}

/*****************************************************************************/

void cm_file_gzopen_read (gzFile *f, const char *filename)
{
   *f = gzopen (filename, "r");
   if (*f == NULL) {
      printf ("Could not open file '%s' for reading.\n", filename);
      exit (1);
   }
}

/*****************************************************************************/

void cm_file_gzclose (gzFile f)
{
   gzclose (f);
}

/*****************************************************************************/

bool cm_class_write (cm_class_srcptr c, cm_param_srcptr param)
   /* Write the class polynomial to the file
      CM_CLASS_DATADIR + "/cp_" + d + "_" + invariant + "_" + paramstr
         + ".dat".
      If an error occurs, the return value is false. */

{
   char filename [400];
   FILE *f;
   int i;

   sprintf (filename, "%s/cp_%"PRIicl"_%c_%s.dat", CM_CLASS_DATADIR, -param->d,
            param->invariant, param->str);

   if (!cm_file_open_write (&f, filename))
      return false;

   fprintf (f, "%"PRIicl"\n", -param->d);
   fprintf (f, "%c\n", param->invariant);
   fprintf (f, "%s\n", param->str);
   fprintf (f, "%i\n", c->classpol->deg);
   for (i = c->classpol->deg; i >= 0; i--) {
      mpz_out_str (f, 10, c->classpol->coeff [i]);
      if (param->field == CM_FIELD_COMPLEX) {
         fprintf (f, " ");
         mpz_out_str (f, 10, c->classpol_c->coeff [i]);
      }
      fprintf (f, "\n");
   }

   cm_file_close (f);

   return true;
}

/*****************************************************************************/

bool cm_class_read (cm_class_ptr c, cm_param_srcptr param)
   /* Read the class polynomial from a file written by cm_class_write.
      If an error occurs, the return value is false. */

{
   char filename [400];
   FILE* f;
   int i;
   char inv;
   char pars [255];
   int_cl_t disc;

   sprintf (filename, "%s/cp_%"PRIicl"_%c_%s.dat", CM_CLASS_DATADIR, -param->d,
            param->invariant, param->str);

   if (!cm_file_open_read (&f, filename))
      return false;

   if (!fscanf (f, "%"SCNicl"\n", &disc))
      return false;
   if (-disc != param->d) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** discriminant %"PRIicl" instead of %"PRIicl"\n",
         -disc, param->d);
      return false;
   }
   if (!fscanf (f, "%c", &inv))
      return false;
   if (inv != param->invariant) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** invariant '%c' instead of '%c'\n",
         inv, param->invariant);
      return false;
   }
   if (!fscanf (f, "%254s", pars))
      return false;
   if (strcmp (pars, param->str)) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** parameter %s instead of %s\n", pars, param->str);
      return false;
   }
   if (!fscanf (f, "%i", &i))
      return false;
   if (i != c->classpol->deg) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** degree %i instead of %i\n", i, c->classpol->deg);
      return false;
   }

   for (i = c->classpol->deg; i >= 0; i--) {
      mpz_inp_str (c->classpol->coeff [i], f, 10);
      if (param->field == CM_FIELD_COMPLEX)
         mpz_inp_str (c->classpol_c->coeff [i], f, 10);
   }

   cm_file_close (f);
   c->computed_classpol = true;

   return true;
}

/*****************************************************************************/

void cm_class_print_pari (FILE* file, cm_class_srcptr c,
   char *fun, char *var)
   /* Print the class polynomial and/or class field tower decomposition
      in c to file using a format understandable by PARI. fun contains
      the function base name ("f" by default), var the base name of the
      variable ("x" by default); the default values are chosen when the
      arguments are NULL. */
{
   char* f = (fun == NULL ? "f" : fun);
   char* x = (var == NULL ? "x" : var);

   if (c->computed_classpol) {
      printf ("%s = ", f);
      if (c->field == CM_FIELD_REAL)
         mpzx_print_pari (file, c->classpol, x);
      else
         mpzxx_print_pari (file, c->classpol, c->classpol_c, x);
      printf ("\n");
   }
   if (c->computed_tower) {
      if (c->field == CM_FIELD_REAL)
         mpzx_tower_print_pari (file, c->tower, f, x);
      else
         mpzxx_tower_print_pari (file, c->tower, c->tower_c, f, x);
   }
}

/*****************************************************************************/

bool cm_file_write_h (const char *tmpdir, const unsigned int *h, unsigned int n)
   /* Write the 2^n class numbers stored in h to a file in tmpdir. */

{
   char *filename;
   FILE *f;
   unsigned int len;
   unsigned long int no;
   bool res;

   len = strlen (tmpdir) + 7;
   filename = (char *) malloc (len * sizeof (char));
   snprintf (filename, len, "%s/h.dat", tmpdir);

   f = fopen (filename, "w");
   res = (f != NULL);

   if (res) {
      no = 1ul << n;
      res = (fwrite (h, sizeof (int), no, f) == no);
      fclose (f);
   }

   free (filename);

   return res;
}

/*****************************************************************************/

bool cm_file_read_h (const char *tmpdir, unsigned int *h, unsigned int n)
   /* Read the 2^n class numbers from a file in tmpdir into h.
      If the file contains fewer class numbers, the reading fails and the
      return value is automatically false. */

{
   char *filename;
   FILE *f;
   unsigned int len;
   unsigned long int no;
   bool res;

   len = strlen (tmpdir) + 7;
   filename = (char *) malloc (len * sizeof (char));
   snprintf (filename, len, "%s/h.dat", tmpdir);

   f = fopen (filename, "r");
   res = (f != NULL);

   if (res) {
      no = 1ul << n;
      res = (fread (h, sizeof (int), no, f) == no);
      fclose (f);
   }

   free (filename);

   return res;
}

/*****************************************************************************/

bool cm_file_write_primorial (const char *tmpdir, mpz_srcptr prim,
   const int i)
   /* Write the prime product prim to a file in tmpdir, assuming it is the
      i-th one in a list. */

{
   char *filename;
   FILE *f;
   unsigned int len;
   bool res;

   len = strlen (tmpdir) + 30;
      /* enough for a 64-bit number i with 20 digits */
   filename = (char *) malloc (len * sizeof (char));
   snprintf (filename, len, "%s/prim_%04i.dat", tmpdir, i);

   f = fopen (filename, "w");
   res = (f != NULL);

   if (res) {
      res = (mpz_out_raw (f, prim) > 0);
      fclose (f);
   }

   free (filename);

   return res;
}

/*****************************************************************************/

bool cm_file_read_primorial (const char *tmpdir, mpz_ptr prim, const int i)
   /* Read the prime product prim from a file in tmpdir, assuming it is the
      i-th one in a list. */

{
   char *filename;
   FILE *f;
   unsigned int len;
   bool res;

   len = strlen (tmpdir) + 30;
      /* enough for a 64-bit number i with 20 digits */
   filename = (char *) malloc (len * sizeof (char));
   snprintf (filename, len, "%s/prim_%04i.dat", tmpdir, i);

   f = fopen (filename, "r");
   res = (f != NULL);

   if (res) {
      res = (mpz_inp_raw (prim, f) > 0);
      fclose (f);
   }

   free (filename);

   return res;
}

/*****************************************************************************/

bool cm_file_write_factor (const char *tmpdir, mpzx_srcptr factor,
   mpzx_srcptr F, mpz_srcptr p)
   /* With factor being a monic factor of the monic polynomial F modulo p,
      write both to a file in tmpdir the name of which is derived from the
      hash of f and p. */

{
   char *filename;
   uint64_t hash;
   FILE *f;
   unsigned int len;
   bool res;

   len = strlen (tmpdir) + 29;
   filename = (char *) malloc (len * sizeof (char));
   hash = mpzx_mod_hash (F, p);
   snprintf (filename, len, "%s/factor_%016"PRIx64".dat", tmpdir, hash);

   f = fopen (filename, "w");
   res = (f != NULL);

   if (res) {
      res &= (mpz_out_str (f, 10, p) != 0);
      res &= (fprintf (f, "\n") != 0);
      res &= (mpzx_out_str (f, 10, F) != 0);
      res &= (fprintf (f, "\n") != 0);
      res &= (mpzx_out_str (f, 10, factor) != 0);
      res &= (fprintf (f, "\n") != 0);
      res &= (fclose (f) == 0);
   }

   free (filename);

   return res;
}

/*****************************************************************************/

static bool write_stat (FILE *f, cm_stat_t stat)
   /* Write the content of stat to the file f. */
{
   int size, i;
   bool ok = true;

   size = sizeof (stat->counter) / sizeof (stat->counter [0]);
   for (i = 0; i < size; i++)
      ok &= (fprintf (f, "%lu ", stat->counter [i]) != 0);
   size = sizeof (stat->timer) / sizeof (stat->timer [0]);
   for (i = 0; i < size - 1; i++)
      ok &= (fprintf (f, "%f %f ",
         stat->timer [i]->elapsed,
         stat->timer [i]->wc_elapsed) != 0);
   ok &= (fprintf (f, "%f %f\n\n",
      stat->timer [size - 1] ->elapsed,
      stat->timer [size - 1] ->wc_elapsed) != 0);

   return ok;
}

/*****************************************************************************/

static bool read_stat (FILE *f, cm_stat_t stat)
   /* Read statistical information from f into stat. */
{
   int size, i;
   bool ok = true;

   size = sizeof (stat->counter) / sizeof (stat->counter [0]);
   for (i = 0; i < size; i++)
      ok &= (fscanf (f, "%lu", stat->counter + i) != 0);
   size = sizeof (stat->timer) / sizeof (stat->timer [0]);
   for (i = 0; i < size - 1; i++)
      ok &= (fscanf (f, "%lf%lf",
         &(stat->timer [i]->elapsed),
         &(stat->timer [i]->wc_elapsed)) != 0);
   ok &= (fscanf (f, "%lf%lf",
      &(stat->timer [size - 1] ->elapsed),
      &(stat->timer [size - 1] ->wc_elapsed)) != 0);

   return ok;
}

/*****************************************************************************/

static bool read_ecpp_cert1_line (FILE *f, mpz_t *line, cm_stat_t stat)
   /* Try to read an additional line of a first step ECPP certificate
      from f and return it in line, which needs to contain four initialised
      entries; the return value indicates the success of the operation.
      Statistical information is also read from the file and returned
      in stat. */
{
   int i;
   bool ok;

   for (i = 0, ok = true; i < 4 && ok; i++)
      ok = (mpz_inp_str (line [i], f, 10) != 0);

   ok &= read_stat (f, stat);

   return ok;
}

/*****************************************************************************/

static int read_ecpp_cert2_line (FILE *f, mpz_t *line, cm_stat_t stat)
   /* Try to read an additional line of a second step ECPP certificate
      from f and return it in line, which needs to contain six initialised
      entries. If the operation was successful, the return value gives the
      number of the entry read; otherwise the return values is -1.
      Statistics information is also read and returned in stat. */
{
   int no, i;
   bool ok;

   ok = (fscanf (f, "%i\n", &no) == 1);

   for (i = 0; i < 6 && ok; i++)
      ok = (mpz_inp_str (line [i], f, 10) != 0);

   ok &= read_stat (f, stat);

   return (ok ? no : -1);
}

/*****************************************************************************/

bool cm_write_ecpp_cert1_line (FILE *f, mpz_t *line, cm_stat_t stat)
   /* Write line, supposed to contain one entry for the first part of ECPP,
      to f; the return value indicates the success of the operation.
      The stat entry is also written so that comprehensive statistics
      of the full run can be collected. */
{
   int i;
   bool ok;

   for (i = 0, ok = true; i < 4 && ok; i++) {
      ok = (mpz_out_str (f, 10, line [i]) != 0);
      ok &= (fprintf (f, "\n") != 0);
   }
   ok &= write_stat (f, stat);

   fflush (f);

   return ok;
}

/*****************************************************************************/

bool cm_write_ecpp_cert2_line (FILE *f, mpz_t *line, int no, cm_stat_t stat)
   /* Write line, supposed to contain the entry number no for the second
      part of ECPP, to f; the return value indicates the success of the
      operation.
      Statistics from stat is also written to the file. */
{
   int i;
   bool ok;

   ok = (fprintf (f, "%i\n", no) != 0);

   for (i = 0; i < 6 && ok; i++) {
      ok = (mpz_out_str (f, 10, line [i]) != 0);
      ok &= (fprintf (f, "\n") != 0);
   }

   ok &= write_stat (f, stat);

   fflush (f);

   return ok;
}

/*****************************************************************************/

mpz_t** cm_file_read_ecpp_cert1 (int *depth, mpz_srcptr p, FILE *f,
   bool debug, cm_stat_t stat)
   /* Try to read a (partial) result of the first ECPP step for p from the
      file f; it is returned via a newly allocated array of size depth,
      and statistics are read from the file and returned in stat.
      The file position indicator is advanced behind the read part, so that
      new entries will be written at the end. */
{
   mpz_t line [4];
   mpz_t **c = NULL;
   int i;

   for (i = 0; i < 4; i++)
      mpz_init (line [i]);
   *depth = 0;
   while (read_ecpp_cert1_line (f, line, stat)) {
      c = (mpz_t **) realloc (c, (*depth + 1) * sizeof (mpz_t *));
      c [*depth] = (mpz_t *) malloc (4 * sizeof (mpz_t));
      for (i = 0; i < 4; i++) {
         c [*depth][i][0] = line [i][0];
         mpz_init (line [i]);
      }
      (*depth)++;
   }
   for (i = 0; i < 4; i++)
      mpz_clear (line [i]);
   if (*depth > 0 && mpz_cmp (p, c [0][0]) != 0) {
      printf ("***** Error: File in cm_file_read_ecpp_cert1 does not "
            "correspond\nto the number to be proved prime.\n");
      exit (1);
   }
   if (debug) {
      printf ("Read %i stage 1 entr", *depth);
      if (*depth == 1)
         printf ("y");
      else
         printf ("ies");
      printf (" from file.\n");
   }

   return c;
}

/*****************************************************************************/

int cm_file_read_ecpp_cert2 (mpz_t **c, mpz_srcptr p, FILE *f, bool debug,
   cm_stat_t stat)
   /* Try to read a (partial) result of the second ECPP step for p from the
      file f; it is returned in c, which needs to be allocated and
      initialised. The return value is the number of read entries.
      Statistical information is also read and returned in stat.
      The file position indicator is advanced behind the read part, so that
      new entries will be written at the end. */
{
   mpz_t *line, *tmp;
   int read, no, i;

   line = (mpz_t *) malloc (6 * sizeof (mpz_t));
   for (i = 0; i < 6; i++)
      mpz_init (line [i]);

   read = 0;
   while ((no = read_ecpp_cert2_line (f, line, stat)) != -1) {
      if (no == 0 && mpz_cmp (p, line [0])) {
         printf ("***** Error: File in cm_file_read_ecpp_cert2 does not "
               "correspond\nto the number to be proved prime.\n");
         exit (1);
      }
      tmp = c [no];
      c [no] = line;
      line = tmp;
      read++;
   }

   for (i = 0; i < 6; i++)
      mpz_clear (line [i]);
   free (line);

   if (debug) {
      printf ("Read %i stage 2 entr", read);
      if (read == 1)
         printf ("y");
      else
         printf ("ies");
      printf (" from file.\n");
   }

   return read;
}

/*****************************************************************************/

void cm_file_write_ecpp_cert_pari (FILE *f, mpz_t **c, int l)
   /* Write the ECPP certificate of lengh l in c to the already opened
      file f, in the format used by PARI. */
{
   int i, j;

   fprintf (f, "[");
   for (i = 0; i < l; i++) {
      fprintf (f, "[");
      for (j = 0; j < 4; j++) {
         mpz_out_str (f, 10, c [i][j]);
         fprintf (f, ", ");
      }
      fprintf (f, "[");
      mpz_out_str (f, 10, c [i][4]);
      fprintf (f, ", ");
      mpz_out_str (f, 10, c [i][5]);
      fprintf (f, "]]");
      if (i != l - 1)
         fprintf (f, ", ");
   }
   fprintf (f, "]\n");
}

/*****************************************************************************/

static void mpz_out_hex (FILE *f, mpz_t z)
   /* Print z to f in hexadecimal format, preceded by 0x and followed
      by a newline character. */
{
   if (mpz_cmp_ui (z, 0) < 0) {
      mpz_neg (z, z);
      fprintf (f, "-0x");
      mpz_out_str (f, 16, z);
      mpz_neg (z, z);
   }
   else {
      fprintf (f, "0x");
      mpz_out_str (f, 16, z);
   }
   fprintf (f, "\n");
}

/*****************************************************************************/

void cm_file_write_ecpp_cert_primo (FILE *f, mpz_t **c, int l)
   /* Write the ECPP certificate of lengh l in c to the already opened
      file f, in the format used by Primo. */
{
   mpz_t tmp, b, j;
   int i;

   mpz_init (tmp);
   mpz_init (b);
   mpz_init (j);

   fprintf (f, "[PRIMO - Primality Certificate]\n"
               "Format=4\n"
               "TestCount=%i\n\n", l);

   fprintf (f, "[Comments]\n"
               "Generated by CM version %s, available under the\n",
               cm_get_version ());
   fprintf (f, "GNU General Public License version 3.0 or later at\n"
               "https://www.multiprecision.org/cm/\n\n");

   fprintf (f, "[Candidate]\nN=");
   mpz_out_hex (f, c [0][0]);
   for (i = 0; i < l; i++) {
      /* The certificate follows the same rules as PARI/GP concerning the
         special cases j=0 or j=1728, but all numbers are output in
         base 10, which seems to work at least with Primo 4.3.3.
         Recall for the following formulae that
         p = c [i][0],
         t = c [i][1],
         co = c [i][2],
         a = c [i][3],
         x = c [i][4],
         y = c [i][5]. */
      fprintf (f, "\n[%i]\nS=", i+1);
      mpz_out_hex (f, c [i][2]);
      fprintf (f, "W=");
      mpz_out_hex (f, c [i][1]);

      /* Compute b. */
      mpz_mul (tmp, c [i][4], c [i][4]);
      mpz_add (tmp, tmp, c [i][3]);
      mpz_mod (tmp, tmp, c [i][0]);
      mpz_mul (tmp, tmp, c [i][4]);
      mpz_mod (tmp, tmp, c [i][0]);
      mpz_mul (b, c [i][5], c [i][5]);
      mpz_sub (b, b, tmp);
      mpz_mod (b, b, c [i][0]);

      /* Compute j. */
      mpz_powm_ui (j, c [i][3], 3, c [i][0]);
      mpz_mul_ui (j, j, 4);
      mpz_mul (tmp, b, b);
      mpz_mul_ui (tmp, tmp, 27);
      mpz_add (tmp, tmp, j);
      mpz_invert (tmp, tmp, c [i][0]);
      mpz_mul (j, j, tmp);
      mpz_mul_ui (j, j, 1728);
      mpz_mod (j, j, c [i][0]);
      /* Centre j around 0. */
      mpz_sub (tmp, j, c [i][0]);
      if (mpz_cmpabs (tmp, j) < 0)
         mpz_set (j, tmp);

      if (!mpz_cmp_ui (j, 0) || !mpz_cmp_ui (j, 1728)) {
         fprintf (f, "A=");
         mpz_out_hex (f, c [i][3]);
         fprintf (f, "B=");
         mpz_out_hex (f, b);
         fprintf (f, "T=");
         mpz_out_hex (f, c [i][4]);
      }
      else {
         fprintf (f, "J=");
         mpz_out_hex (f, j);
         /* Compute T = 2 * (1728 - j) * a * x / (3 * b). */
         mpz_ui_sub (tmp, 1728, j);
         mpz_mul_ui (tmp, tmp, 2);
         mpz_mul (tmp, tmp, c [i][3]);
         mpz_mod (tmp, tmp, c [i][0]);
         mpz_mul (tmp, tmp, c [i][4]);
         mpz_mod (tmp, tmp, c [i][0]);
         mpz_mul_ui (b, b, 3);
         mpz_invert (b, b, c [i][0]);
         mpz_mul (tmp, tmp, b);
         mpz_mod (tmp, tmp, c [i][0]);
         fprintf (f, "T=");
         mpz_out_hex (f, tmp);
      }
   }

   mpz_clear (tmp);
   mpz_clear (b);
   mpz_clear (j);
}

/*****************************************************************************/
