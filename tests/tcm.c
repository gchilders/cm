/*

tcm.c - tests for cm

Copyright (C) 2009, 2010, 2011, 2012, 2016, 2021 Andreas Enge

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

#include <string.h>
#include "cm.h"

/*****************************************************************************/

static void test_curve (int_cl_t d, char invariant, bool verbose)
{
   cm_param_t param;
   cm_class_t c;
   mpz_t a, b, x, y, p, n, l, co;
   cm_timer clock;

   cm_timer_start (clock);

   if (!cm_param_init (param, d, invariant, -1, verbose))
      exit (1);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);
   mpz_init (p);
   mpz_init (n);
   mpz_init (l);
   mpz_init (co);

   cm_curve_crypto_param (p, n, l, co, d, 200, verbose);

   /* First test with class polynomials. */
   printf ("d = %"PRIicl", inv = %c; ", d, invariant);
   printf ("class polynomial: ");
   fflush (stdout);
   cm_class_init (c, param, verbose);
   cm_class_compute (c, param, true, false, verbose);
   cm_curve_and_point (a, b, x, y, param, c, p, l, co, CM_MODPOLDIR,
      verbose);
      /* CM_MODPOLDIR is a preprocessor variable defined by the -D
         parameter of gcc */
   cm_class_clear (c);
   printf ("ok; ");

   /* Then test with a class field tower. */
   printf ("class field tower: ");
   fflush (stdout);
   cm_class_init (c, param, verbose);
   cm_class_compute (c, param, false, true, verbose);
   cm_curve_and_point (a, b, x, y, param, c, p, l, co, CM_MODPOLDIR,
      verbose);
   cm_class_clear (c);
   printf ("ok\n");

   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (x);
   mpz_clear (y);
   mpz_clear (p);
   mpz_clear (n);
   mpz_clear (l);
   mpz_clear (co);

   cm_timer_stop (clock);
   if (verbose)
      printf ("\n--- Elapsed time: %.1f\n", cm_timer_get (clock));
}


/*****************************************************************************/

static void small_test (void)
{
   test_curve ((int_cl_t) (-108715), CM_INVARIANT_DOUBLEETA, false);
}

/*****************************************************************************/

static void big_test (void)
{
   /* Weber: d divisible by 4, not by 32, not by 3 */
   test_curve ((int_cl_t) (-108740), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108712), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108716), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108752), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108724), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108728), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108764), CM_INVARIANT_WEBER, false);
   /* Weber: d divisible by 12, not by 32 */
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108720), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108732), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108744), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108756), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108780), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108792), CM_INVARIANT_WEBER, false);

   /* d==-16, corresponds to fixed bug */
   test_curve ((int_cl_t) (-16), CM_INVARIANT_J, false);

   test_curve ((int_cl_t) (-108708), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108703), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_GAMMA2, false);
   test_curve ((int_cl_t) (-108703), CM_INVARIANT_GAMMA2, false);
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_GAMMA3, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_ATKIN, false); /* p=71 */
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_ATKIN, false); /* p=47 */
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_ATKIN, false); /* p=59 */
   test_curve ((int_cl_t) (-58767),  CM_INVARIANT_ATKIN, false); /* p=131 */
   test_curve ((int_cl_t) (-299),    CM_INVARIANT_SIMPLEETA, false);
   test_curve ((int_cl_t) (-105131), CM_INVARIANT_MULTIETA, false); /* N=3*5*7 */

   test_curve ((int_cl_t) (-1032), CM_INVARIANT_DOUBLEETA, false);
      /* N=3*7, s=e=2 */
}

/*****************************************************************************/

int main (void)
{
   cm_pari_init ();

   small_test ();
   big_test ();
   
   cm_pari_clear ();

   return 0;
}
