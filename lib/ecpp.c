/*

ecpp.c - code for computing ECPP certificates

Copyright (C) 2021 Andreas Enge

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

#include "cm-impl.h"


/*****************************************************************************/

void cm_ecpp (mpz_srcptr N, const char* modpoldir, bool tower, bool print,
   bool verbose)
   /* Assuming that N is a (probable) prime, compute an ECPP certificate.
      modpoldir gives the directory where modular polynomials are stored;
      it is passed through to the function computing a curve from a root
      of the class polynomial.
      tower indicates whether a class field tower decomposition is used
      instead of only the class polynomial.
      print indicates whether the result is printed.
      verbose indicates whether intermediate computations output progress
      information. */
{
   int depth;
   mpz_t **cert;
   mpz_srcptr p, n, l;
   int_cl_t d;
   mpz_t t, co, a, b, x, y;
   cm_param_t param;
   cm_class_t c;
   int i, j;
   cm_timer clock, clock1, clock2;

   mpz_init (t);
   mpz_init (co);
   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);

   cm_timer_start (clock);
   cert = cm_pari_ecpp1 (&depth, N);
   cm_timer_stop (clock);
   if (verbose)
      printf ("--- Time for first ECPP step:  %.1f\n", cm_timer_get (clock));

   cm_timer_start (clock);
   cm_timer_reset (clock1);
   cm_timer_reset (clock2);
   if (print)
      printf ("c = [");
   for (i = 0; i < depth; i++) {
      p = cert [i][0];
      d = mpz_get_si (cert [i][1]);
      n = cert [i][2];
      l = cert [i][3];

      mpz_add_ui (t, p, 1);
      mpz_sub (t, t, n);
      mpz_divexact (co, n, l);

      cm_timer_continue (clock1);
      if (d % 3 != 0)
         cm_param_init (param, d, CM_INVARIANT_GAMMA2, false);
      else if (d % 2 != 0)
         cm_param_init (param, d, CM_INVARIANT_GAMMA3, false);
      else
         cm_param_init (param, d, CM_INVARIANT_J, false);
      cm_class_init (c, param, false);
      cm_class_compute (c, param, !tower, tower, false);
      cm_timer_stop (clock1);
      cm_timer_continue (clock2);
      cm_curve_and_point (a, b, x, y, param, c, p, l, co,
         modpoldir, false);
      cm_timer_stop (clock2);
      cm_class_clear (c);

      if (print) {
         printf ("[");
         mpz_out_str (stdout, 10, p);
         printf (", ");
         mpz_out_str (stdout, 10, t);
         printf (", ");
         mpz_out_str (stdout, 10, co);
         printf (", ");
         mpz_out_str (stdout, 10, a);
         printf (", [");
         mpz_out_str (stdout, 10, x);
         printf (", ");
         mpz_out_str (stdout, 10, y);
         printf ("]]");
         if (i != depth - 1)
            printf (", ");
      }
   }
   if (print)
      printf ("];\n");

   mpz_clear (t);
   mpz_clear (co);
   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (x);
   mpz_clear (y);
   for (i = 0; i < depth; i++) {
      for (j = 0; j < 4; j++)
         mpz_clear (cert [i][j]);
      free (cert [i]);
   }
   free (cert);
   cm_timer_stop (clock);
   if (verbose) {
      printf ("--- Time for CM:               %.1f\n", cm_timer_get (clock1));
      printf ("--- Time for curves:           %.1f\n", cm_timer_get (clock2));
      printf ("--- Time for second ECPP step: %.1f\n", cm_timer_get (clock));
   }
}

/*****************************************************************************/
/*****************************************************************************/
