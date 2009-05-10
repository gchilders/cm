#include "cm.h"

/*****************************************************************************/

static bool evaluate_parameters (int argc, char* argv [], int_cl_t *d,
      char *invariant)

   /* The function determines the parameter values and consequently sets the */
   /* discriminant d, the precision and the type of class invariant. If an   */
   /* error occurs which forces the programme to break the return value is   */
   /* false.                                                                 */

{  int  index = 1;
      /* points to the currently considered entry of the parameter list */
   bool ok = true;
   char *invariant_string = NULL;

   *d = 0;
   *invariant = CM_INVARIANT_NONE;

   while (index < argc && ok)
   {  /* analyse entry "index" of argv */

      if (argv [index] [0] != '-')
      {  printf ("Options must begin with a '-'.\n");
         ok = false;
      }
      else /* entry is an option, check for type */

      if (strlen (argv [index]) >= 2 && argv [index] [1] == 'd')
      {  if (strlen (argv [index]) == 2)
         {  printf ("You specified the option '-d' without any integer following; it should be\n");
            printf ("followed by the absolute value of the discriminant.\n");
            ok = false;
         }
         if (*d != 0)
         {  printf ("You specified both the options '-d%"PRIicl"' and '%s';", -(*d), argv [index]);
            printf ("please decide for one of them.\n");
            ok = false;
         }
         else
         {
            *d = - atoll (argv [index] + 2);
            index ++;
         }
      }

      else if (strlen (argv [index]) >= 2 && argv [index] [1] == 'i')
      {
         if (strlen (argv [index]) == 2)
         {  printf ("You specified the option '-i' without anything following; it should be\n");
            printf ("followed by 'j', 'gamma2', 'gamma3', 'weber', 'doubleeta'\n");
            printf ("or 'simpleeta', depending on which type of class invariant you would\n");
            printf ("like to use for the computations.\n");
            ok = false;
         }
         else if (*invariant != CM_INVARIANT_NONE)
         {  printf ("You specified both the options '-i%s' and '%s';", invariant_string, argv [index]);
            printf ("please decide for one of them.\n");
            ok = false;
         }
         else
         {
            invariant_string = argv [index] + 2;
            if      (!strcmp (invariant_string, "j"))
               *invariant = CM_INVARIANT_J;
            else if (!strcmp (invariant_string, "gamma2"))
               *invariant = CM_INVARIANT_GAMMA2;
            else if (!strcmp (invariant_string, "gamma3"))
               *invariant = CM_INVARIANT_GAMMA3;
            else if (!strcmp (invariant_string, "weber"))
               *invariant = CM_INVARIANT_WEBER;
            else if (!strcmp (invariant_string, "doubleeta"))
               *invariant = CM_INVARIANT_DOUBLEETA;
            else if (!strcmp (invariant_string, "simpleeta"))
               *invariant = CM_INVARIANT_SIMPLEETA;
            else
            {
               printf ("You specified the option '-i' follow by '%s',", invariant_string);
               printf ("which is not a recognised\n");
               printf ("class invariant. It should be followed by 'j', 'gamma2', 'gamma3', 'weber',\n");
               printf ("'doubleeta' or 'simpleeta'\n");
               ok = false;
            }
            index++;
         }
      }

      else
      {  printf ("You specified the option '%s' ", argv [index]);
         printf ("which does not exist. You may use the options\n");
         printf ("'-d' followed by the absolute value of the discriminant or\n");
         printf ("'-i' followed by 'j', 'gamma2', 'gamma3', 'weber', 'doubleeta'\n");
         printf ("     or 'simpleeta', depending on which type of class invariant you\n");
         printf ("     would like to use for the computations.\n");
         ok = false;
      }

   }

   if (*d >= 0 || (*d % 4 != 0 && (*d - 1) % 4 != 0)) {
      printf ("d = %"PRIicl" is not a quadratic discriminant\n", *d);
      ok = false;
   }

   return ok;

}

/*****************************************************************************/

int main (int argc, char* argv [])
{
   char     invariant;
   int_cl_t d;
   cm_timer clock;

   cm_timer_start (clock);

   if (VERBOSE)
      printf ("_______________________________________________________________"
              "________________\n\n");
   if (!evaluate_parameters (argc, argv, &d, &invariant))
      exit (1);

   if (invariant == CM_INVARIANT_NONE)
      invariant = CM_INVARIANT_DOUBLEETA;
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
