#include "cm_class.h"
#include "string.h"

/*****************************************************************************/

void test_curve (int_cl_t d, char invariant, bool verbose) {
   cm_timer clock;

   cm_timer_start (clock);

   printf ("d = %"PRIicl", inv = %c ", d, invariant);
   fflush (stdout);
   cm_curve_compute_curve (d, invariant, CM_MODPOLDIR, verbose);
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
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_DOUBLEETA, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_J, false);
   test_curve ((int_cl_t) (-108708), CM_INVARIANT_WEBER, false);
   test_curve ((int_cl_t) (-108707), CM_INVARIANT_GAMMA2, false);
   test_curve ((int_cl_t) (-108711), CM_INVARIANT_GAMMA3, false);
   test_curve ((int_cl_t) (-299), CM_INVARIANT_SIMPLEETA, false);
//   test_curve ((int_cl_t) (-6961631), CM_INVARIANT_DOUBLEETA, true);

   return 0;
}
