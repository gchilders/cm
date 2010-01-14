/*

tcm.c - tests for cm

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

#include "cm_class.h"
#include "string.h"

/*****************************************************************************/

void test_curve (int_cl_t d, char invariant, bool verbose) {
   cm_timer clock;

   cm_timer_start (clock);

   printf ("d = %"PRIicl", inv = %c ", d, invariant);
   fflush (stdout);
   cm_curve_compute_curve (d, invariant, 200, CM_MODPOLDIR, false, verbose);
      /* CM_MODPOLDIR is a preprocessor variable defined by the -D
         parameter of gcc */
   printf (" ok\n");

   cm_timer_stop (clock);
   if (verbose)
      printf ("\n--- Elapsed time: %.1f\n", cm_timer_get (clock));
}


/*****************************************************************************/

int main ()
{
   test_curve ((int_cl_t) (-108715), CM_INVARIANT_DOUBLEETA, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_GAMMA2, false);
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_GAMMA3, false);
   test_curve ((int_cl_t) (-299), CM_INVARIANT_SIMPLEETA, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_ATKIN, false); /* p=71 */
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_ATKIN, false); /* p=47 */
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_ATKIN, false); /* p=59 */

   return 0;
}
