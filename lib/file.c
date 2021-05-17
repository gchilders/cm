/*

file.c - code for handling (gzipped) files

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

/*****************************************************************************/

bool cm_file_open_write (FILE **f, char *filename)
{
   *f = fopen (filename, "w");
   if (*f == NULL) {
      printf ("Could not open file '%s' for writing.\n", filename);
      return false;
   }
   else {
      printf ("Writing to '%s'.\n", filename);
      return true;
   }
}

/*****************************************************************************/

bool cm_file_open_read (FILE **f, char *filename)
{
   *f = fopen (filename, "r");
   if (*f == NULL) {
      printf ("Could not open file '%s' for reading.\n", filename);
      return false;
   }
   else {
      printf ("Reading from '%s'.\n", filename);
      return true;
   }
}

/*****************************************************************************/

void cm_file_close (FILE *f)
{
   fclose (f);
}

/*****************************************************************************/

void cm_file_gzopen_write (gzFile *f, char *filename)
{
   *f = gzopen (filename, "w9");
   if (*f == NULL) {
      printf ("Could not open file '%s' for writing.\n", filename);
      exit (1);
   }
}

/*****************************************************************************/

void cm_file_gzopen_read (gzFile *f, char *filename)
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
   char *fun, char *fun_c, char *var)
   /* Print the class polynomial and/or class field tower decomposition in c
      to file using a format understandable by PARI. fun and fun_c contain
      the function names for the real and, if applicable, the complex parts
      of the tower ("f" and "g" by default), var the base name of the
      variable ("x" by default); the default values are chosen when the
      arguments are NULL. */
{
   char* f = (fun == NULL ? "f" : fun);
   char* g = (fun_c == NULL ? "g" :fun_c);
   char* x = (var == NULL ? "x" : var);

   if (c->computed_classpol) {
      printf ("%s = ", f);
      if (c->field == CM_FIELD_REAL) {
         mpzx_print_pari (file, c->classpol, x);
         printf (")\n");
      }
      else {
         printf ("(");
         mpzx_print_pari (file, c->classpol, x);
         printf (")+o*(");
         mpzx_print_pari (file, c->classpol_c, x);
         printf (")\n");
      }
   }
   if (c->computed_tower) {
      if (c->field == CM_FIELD_REAL)
         mpzx_tower_print_pari (file, c->tower, f, x);
      else {
         mpzx_tower_print_pari (file, c->tower, f, x);
         mpzx_tower_print_pari (file, c->tower_c, g, x);
      }
   }
}

/*****************************************************************************/
