/*

cm.c - executable computing a cryptographically suitable cm curve

Copyright (C) 2009, 2010, 2021 Andreas Enge

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

int main (int argc, char* argv [])
{
   int_cl_t d;
   char     invariant;
   bool     verbose;
   cm_param_t param;
   cm_class_t c;
   mpz_t a, b, x, y, p, n, l, co;
   cm_timer clock;

   cm_timer_start (clock);

   if (!evaluate_parameters (argc, argv, &d, &invariant, &verbose))
      exit (1);
   if (invariant == CM_INVARIANT_NONE)
      invariant = CM_INVARIANT_J;

   if (!cm_param_init (param, d, invariant,
      -1, CM_SUBFIELD_OPTIMAL, verbose))
      exit (1);
   cm_pari_init ();
   cm_class_init (c, param, verbose);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);
   mpz_init (p);
   mpz_init (n);
   mpz_init (l);
   mpz_init (co);

   cm_class_compute (c,
      param,
      false /* classpol */,
      true /* tower */,
      verbose);
   cm_curve_crypto_param (p, n, l, co, d, 256, verbose);
   cm_curve_and_point (a, b, x, y, param, c, p, l, co, CM_MODPOLDIR,
      verbose);
      /* CM_MODPOLDIR is a preprocessor variable defined by the -D
         parameter of gcc */

   cm_class_clear (c);
   cm_pari_clear ();
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

   return 0;
}
