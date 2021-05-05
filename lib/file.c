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
