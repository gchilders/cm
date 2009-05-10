#include "cm_class.h"
#include "string.h"

#define VERBOSE false

/*****************************************************************************/

int main ()
{
   char     invariant;
   int_cl_t d;
   cm_timer clock;

   cm_timer_start (clock);

   invariant = CM_INVARIANT_DOUBLEETA;
   d = (int_cl_t) (-108708);

   if (VERBOSE) {
      printf ("Using ");
      switch (invariant)
      {
      case CM_INVARIANT_J:
         printf ("j");
         break;
      case CM_INVARIANT_GAMMA2:
         printf ("gamma2");
         break;
      case CM_INVARIANT_GAMMA3:
         printf ("gamma3");
         break;
      case CM_INVARIANT_WEBER:
         printf ("Weber's function");
         break;
      case CM_INVARIANT_DOUBLEETA:
         printf ("a double eta quotient");
         break;
      case CM_INVARIANT_SIMPLEETA:
         printf ("a simple eta quotient");
         break;
      }
      printf (" as class invariant.\n");

      printf ("d = %"PRIicl"\n", d);
      if (invariant == CM_INVARIANT_WEBER)
      {
         if (d % 3 == 0)
            printf ("d is divisible by 3, and ");
         if ((d - 1) % 8 == 0)
            printf ("d is congruent to 1 modulo 8.\n");
         else if ((d - 5) % 8 == 0)
            printf ("d is congruent to 5 modulo 8.\n");
         else
            printf ("d/4 is congruent to %"PRIicl" modulo 8\n", (d/4) % 8 + 8);
      }
      printf ("\n");
   }

   cm_curve_compute_curve (d, invariant, CM_MODPOLDIR, VERBOSE);
      /* CM_MODPOLDIR is a preprocessor variable defined by the -D
         parameter of gcc */

   cm_timer_stop (clock);
   if (VERBOSE)
      printf ("\n--- Elapsed time: %.1f\n", cm_timer_get (clock));

   return 0;
}
