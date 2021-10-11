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

static void print_d_options (void);
static void print_i_options (void);
static void print_n_options (void);
static void print_v_options (void);
static void print_o_options (void);
static void print_g_options (void);
static void print_help (void);
static void print_help_ecpp (void);
static void print_libraries (void);

/*****************************************************************************/

static void print_d_options (void)
{
   printf ("-d followed by the absolute value of the discriminant "
      "is a required\n"
      "   parameter.\n");
}

/*****************************************************************************/

static void print_i_options (void)
{
   printf ("-i should be followed by one of the following selections:\n"
           "   'j', 'gamma2', 'gamma3', 'weber', 'doubleeta', "
           "'simpleeta',\n"
           "   'multieta' or 'atkin'\n");
}

/*****************************************************************************/

static void print_n_options (void)
{
   printf ("-n followed by a positive number to be proved prime is a "
      "required parameter.\n");
}

/*****************************************************************************/

static void print_v_options (void)
{
   printf ("-v enables verbose output.\n");
}

/*****************************************************************************/

static void print_o_options (void)
{
   printf ("-o enables output of the certificate.\n");
}

/*****************************************************************************/

static void print_g_options (void)
{
   printf ("-g enables debug output.\n");
}

/*****************************************************************************/

static void print_help (void)
{
   printf ("The following options are recognised: "
      "'-d', '-i', '-v', '-h'.\n"
      "-h prints this help.\n");
   print_d_options ();
   print_i_options ();
   print_v_options ();
}

/*****************************************************************************/

static void print_help_ecpp (void)
{
   printf ("The following options are recognised: "
      "'-n', '-o', '-v', '-g', '-h'.\n"
      "-h prints this help.\n");
   print_n_options ();
   print_o_options ();
   print_v_options ();
   print_g_options ();
}

/*****************************************************************************/

static void print_libraries (void)
{
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

/*****************************************************************************/

void evaluate_parameters (int argc, char* argv [], int_cl_t *d,
   char *invariant, bool *verbose)
   /* The function determines the parameter values and consequently sets
      the discriminant d, the type of class invariant and the verbose
      parameter. */
{
   int opt;
   char *invariant_string = NULL;

   *d = 0;
   *invariant = CM_INVARIANT_NONE;
   *verbose = false;

   while ((opt = getopt (argc, argv, "hd:i:v")) != -1) {
      switch (opt) {
         case 'v':
            *verbose = true;
            break;
         case 'd':
            *d = - atoll (optarg);
            if (*d >= 0 || (*d % 4 != 0 && (*d - 1) % 4 != 0)) {
               printf ("d = %"PRIicl" is not a negative quadratic "
                  "discriminant\n", *d);
               exit (1);
            }
            break;
         case 'i':
            invariant_string = optarg;
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
               printf ("You specified the option '-i' followed by '%s'.\n",
                  invariant_string);
               print_i_options ();
               exit (1);
            }
            break;
         case 'h':
            print_help ();
            exit (0);
         case '?':
            if (optopt == 'i')
               print_i_options ();
            else if (optopt == 'd')
               print_d_options ();
            else if (isprint (optopt)) {
               printf ("Unknown option '-%c'.\n", optopt);
               print_help ();
            }
            else {
               printf ("Unknown option with character code %i.\n", optopt);
               print_help ();
            }
            exit (1);
         default:
            /* Should not occur. */
            exit (1);
      }
   }

   if (*d == 0) {
      print_d_options ();
      exit (1);
   }

   if (*verbose)
      print_libraries ();
}

/*****************************************************************************/

void evaluate_parameters_ecpp (int argc, char* argv [], mpz_ptr n,
   bool *output, bool *verbose, bool *debug)
   /* Since ECPP requires different parameter types, the easiest solution
      appears to be a separate function, albeit with a lot of copy and
      paste. */
{
   int opt;

   mpz_set_ui (n, 0ul);
   *output = false;
   *verbose = false;
   *debug = false;

   while ((opt = getopt (argc, argv, "hn:ogv")) != -1) {
      switch (opt) {
         case 'v':
            *verbose = true;
            break;
         case 'o':
            *output = true;
            break;
         case 'g':
            *verbose = true;
            *debug = true;
            break;
         case 'n':
            mpz_set_str (n, optarg, 10);
            break;
         case 'h':
            print_help_ecpp ();
            exit (0);
         case '?':
            if (optopt == 'n')
               print_n_options ();
            else if (isprint (optopt)) {
               printf ("Unknown option '-%c'.\n", optopt);
               print_help_ecpp ();
            }
            else {
               printf ("Unknown option with character code %i.\n", optopt);
               print_help_ecpp ();
            }
            exit (1);
         default:
            /* Should not occur. */
            exit (1);
      }
   }

   if (mpz_cmp_ui (n, 0ul) <= 0) {
      print_n_options ();
      exit (1);
   }

   if (*verbose)
      print_libraries ();
}

/*****************************************************************************/
