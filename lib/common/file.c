#include "cm_common-impl.h"

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

void cm_file_gzopen_write (FILE **f, char *filename)
{
   *f = gzopen (filename, "w9");
   if (*f == NULL) {
      printf ("Could not open file '%s' for writing.\n", filename);
      exit (1);
   }
}

/*****************************************************************************/

void cm_file_gzopen_read (FILE **f, char *filename)
{
   *f = gzopen (filename, "r");
   if (*f == NULL) {
      printf ("Could not open file '%s' for reading.\n", filename);
      exit (1);
   }
}

/*****************************************************************************/

void cm_file_gzclose (FILE *f)
{
   gzclose (f);
}

/*****************************************************************************/
