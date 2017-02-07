/*

params.c - command line parameter evaluation

Copyright (C) 2009, 2010 Andreas Enge

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

#include "params.h"

/*****************************************************************************/

bool evaluate_parameters (int argc, char* argv [], int_cl_t *d,
   char *invariant, bool *verbose)

   /* The function determines the parameter values and consequently sets the */
   /* discriminant d, the type of class invariant and the verbose parameter. */
   /* If an error occurs which forces the program to break the return value  */
   /* is false.                                                              */

{
   int  index = 1;
      /* points to the currently considered entry of the parameter list */
   bool ok = true, paramd = false;
   char *invariant_string = NULL;

   *d = 0;
   *invariant = CM_INVARIANT_NONE;
   *verbose = false;

   while (index < argc && ok) {
      /* analyse entry "index" of argv */

      if (argv [index] [0] != '-') {
         printf ("Options must begin with a '-'.\n");
         ok = false;
      }
      else /* entry is an option, check for type */

      if (strlen (argv [index]) == 2 && argv [index] [1] == 'v') {
         *verbose = true;
         index++;
      }

      else if (strlen (argv [index]) >= 2 && argv [index] [1] == 'd') {
         if (strlen (argv [index]) == 2) {
            printf ("You specified the option '-d' without any integer following; it should be\n");
            printf ("followed by the absolute value of the discriminant.\n");
            ok = false;
         }
         if (*d != 0) {
            printf ("You specified both the options '-d%"PRIicl"' and '%s';", -(*d), argv [index]);
            printf ("please decide for one of them.\n");
            ok = false;
         }
         else {
            *d = - atoll (argv [index] + 2);
            paramd = true;
            index++;
         }
      }

      else if (strlen (argv [index]) >= 2 && argv [index] [1] == 'i') {
         if (strlen (argv [index]) == 2) {
            printf ("You specified the option '-i' without anything following; it should be\n");
            printf ("followed by 'j', 'gamma2', 'gamma3', 'weber', 'doubleeta', 'simpleeta',\n");
            printf ("'multieta' or 'atkin', depending on which type of class invariant you\n");
            printf ("would like to use for the computations.\n");
            ok = false;
         }
         else if (*invariant != CM_INVARIANT_NONE) {
            printf ("You specified both the options '-i%s' and '%s';", invariant_string, argv [index]);
            printf ("please decide for one of them.\n");
            ok = false;
         }
         else {
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
            else if (!strcmp (invariant_string, "multieta"))
               *invariant = CM_INVARIANT_MULTIETA;
            else if (!strcmp (invariant_string, "atkin"))
               *invariant = CM_INVARIANT_ATKIN;
            else {
               printf ("You specified the option '-i' follow by '%s', ", invariant_string);
               printf ("which is not a recognised\n");
               printf ("class invariant. It should be followed by 'j', 'gamma2', 'gamma3', 'weber',\n");
               printf ("'doubleeta', 'simpleeta', 'multieta' or 'atkin'\n");
               ok = false;
            }
            index++;
         }
      }

      else {
         printf ("You specified the option '%s' ", argv [index]);
         printf ("which does not exist. You may use the options\n");
         printf ("'-d' followed by the absolute value of the discriminant\n");
         printf ("'-i' followed by 'j', 'gamma2', 'gamma3', 'weber', 'doubleeta',\n");
         printf ("     'simpleeta', 'multieta' or 'atkin', depending on which type\n");
         printf ("     of class invariant you would like to use for the computations\n");
         printf ("'-v' to enable verbose output\n");
         ok = false;
      }

   }

   if (!paramd) {
      printf ("Please specify '-d', followed by the absolute value of the discriminant.\n");
      ok = false;
   }
   else if (*d >= 0 || (*d % 4 != 0 && (*d - 1) % 4 != 0)) {
      printf ("d = %"PRIicl" is not a quadratic discriminant\n", *d);
      ok = false;
   }

   if (*verbose) {
      GEN v;
      printf ("GMP: include %d.%d.%d, lib %s\n",
            __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL,
            gmp_version);
      printf ("MPFR: include %s, lib %s\n",
            MPFR_VERSION_STRING,
            mpfr_get_version ());
      printf ("MPC: include %s, lib %s\n", MPC_VERSION_STRING,
            mpc_get_version ());
      printf ("MPFRCX: include %s, lib %s\n", MPFRCX_VERSION_STRING,
            mpfrcx_get_version ());
      pari_init (100000, 0);
      v = pari_version ();
      printf ("PARI: include %i.%li.%li, lib %li.%li.%li\n",
            PARI_VERSION_CODE >> 16, (PARI_VERSION_CODE >> 8) & 255ul,
            PARI_VERSION_CODE & 255ul,
            itos (gel (v, 1)), itos (gel (v, 2)), itos (gel (v, 3)));
      pari_close ();
   }

   return ok;
}

/*****************************************************************************/
