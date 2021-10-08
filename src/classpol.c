/*

classpol.c - executable computing a ring class polynomial

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
   int_cl_t   d;
   char       invariant;
   bool       verbose;
   cm_class_t c;
   cm_timer clock;
   cm_param_t param;

   cm_timer_start (clock);

   evaluate_parameters (argc, argv, &d, &invariant, &verbose);
   if (invariant == CM_INVARIANT_NONE)
      invariant = CM_INVARIANT_J;

   if (!cm_param_init (param, d, invariant,
      0, CM_SUBFIELD_NEVER, verbose))
      exit (1);
   cm_pari_init ();
   cm_class_init (c, param, verbose);
   cm_class_compute (c,
      param,
      true /* classpol */,
      false /* tower */,
      verbose);
   cm_class_print_pari (stdout, c, NULL, NULL);
   cm_class_clear (c);
   cm_pari_clear ();

   cm_timer_stop (clock);
   if (verbose)
      printf ("\n--- Elapsed time: %.1f\n", cm_timer_get (clock));

   return 0;
}
