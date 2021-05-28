/*

param.c - code for handling CM parameters

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

#include "cm-impl.h"

static bool doubleeta_compute_parameter (cm_param_ptr param, int_cl_t d);

/*****************************************************************************/

bool cm_param_init (cm_param_ptr param, int_cl_t d, char invariant,
   bool verbose)
   /* Test whether the discriminant is suited for the chosen invariant and
      in this case compute and store the parameters in param and return
      true; otherwise return false. */

{
   int i;
   char* pointer;

   if (d >= 0) {
      printf ("\n*** The discriminant must be negative.\n");
      exit (1);
   }
   else if (d % 4 != 0 && (d - 1) % 4 != 0) {
      printf ("\n*** %"PRIicl" is not a quadratic discriminant.\n", d);
      exit (1);
   }

   param->d = d;
   param->invariant = invariant;

   switch (invariant) {
      case CM_INVARIANT_J:
         param->p [0] = 1;
         param->p [1] = 0;
         param->s = 1;
         param->e = 1;
         break;
      case CM_INVARIANT_GAMMA2:
         if (d % 3 == 0) {
            if (verbose) {
               printf ("\n*** %"PRIicl" is divisible by 3, so that gamma2 ",
                        d);
               printf ("cannot be used.\n");
            }
            return false;
         }
         param->p [0] = 1;
         param->p [1] = 0;
         param->s = 1;
         param->e = 1;
         break;
      case CM_INVARIANT_GAMMA3:
         if (d % 2 == 0) {
            if (verbose) {
               printf ("\n*** %"PRIicl" is divisible by 4, so that gamma3 ",
                        d);
               printf ("cannot be used.\n");
            }
            return false;
         }
         param->p [0] = 1;
         param->p [1] = 0;
         param->s = 1;
         param->e = 1;
         break;
      case CM_INVARIANT_WEBER:
         if (d % 32 == 0) {
            if (verbose) {
               printf ("\n*** %"PRIicl" is divisible by 32, so that the Weber ",
                        d);
               printf ("functions cannot be used.\n");
            }
            return false;
         }
         else if (cm_classgroup_mod (d, 8) == 5) {
            if (verbose)
               printf ("\n*** %"PRIicl" is 5 mod 8; the Weber function "
                  "cannot be used directly, but you\nmay compute the ring "
                  "class field for the discriminant %"PRIicl" and apply\n"
                  "a 2-isogeny when computing the elliptic curve.\n",
                  d, 4*d);
            return false;
         }
         else if (cm_classgroup_mod (d, 8) == 1) {
            if (verbose)
               printf ("\n*** %"PRIicl" is 1 mod 8; the Weber function "
                  "cannot be used directly, but you\nmay compute the ring "
                  "class field for the discriminant %"PRIicl", which is\n"
                  "the same as the Hilbert class field\n",
                  d, 4*d);
            return false;
         }

         /* Let m = -disc/4 and p [0] = m % 8. */
         param->p [0] = ((-d) / 4) % 8;
         param->p [1] = 0;
         param->s = 24;
         param->e = 1;
         if (param->p [0] == 1 || param->p [0] == 2 || param->p [0] == 6)
            param->e *= 2;
         else if (param->p [0] == 4 || param->p [0] == 5)
            param->e *= 4;
         if (d % 3 == 0)
            param->e *= 3;
         break;
      case CM_INVARIANT_ATKIN:
         if (cm_nt_kronecker (d, (int_cl_t) 71) != -1)
            /* factor 36, T_5 + T_29 + 1 */
            param->p [0] = 71;
         else if (cm_nt_kronecker (d, (int_cl_t) 131) != -1)
            /* factor 33, T_61 + 1 */
            param->p [0] = 131;
         else if (cm_nt_kronecker (d, (int_cl_t) 59) != -1)
            /* factor 30, T_5 + T_29 */
            param->p [0] = 59;
         else if (cm_nt_kronecker (d, (int_cl_t) 47) != -1)
            /* factor 24, -T_17 */
            param->p [0] = 47;
         else {
            if (verbose) {
               printf ("\n*** 47, 59, 71 and 131 are inert for %"PRIicl, d);
               printf (", so that atkin cannot be used.\n");
            }
            return false;
         }
         param->p [1] = 0;
         param->s = 1;
         param->e = 1;
         break;
      case CM_INVARIANT_MULTIETA:
         param->p [0] = 3;
         param->p [1] = 5;
         param->p [2] = 7;
         param->p [3] = 0;
         param->s = 1;
         param->e = 1;
         break;
      case CM_INVARIANT_SIMPLEETA:
         param->p [0] = 3;
         if (cm_nt_kronecker (d, (int_cl_t) (param->p [0])) == -1) {
            if (verbose)
               printf ("*** Unsuited discriminant\n\n");
            return false;
         }

         param->p [1] = 0;
         param->s = 12;
         param->e = 12;
         break;
      case CM_INVARIANT_DOUBLEETA:
         if (!doubleeta_compute_parameter (param, d))
            return false;
         break;
      default: /* should not occur */
         printf ("class_compute_parameter called for "
                  "unknown class invariant '%c'\n", invariant);
         exit (1);
   }

   /* create parameter string */
   pointer = param->str;
   i = 0;
   do {
      pointer += sprintf (pointer, "%i_", param->p [i]);
      i++;
   } while (param->p [i] != 0);
   sprintf (pointer, "%i_%i", param->e, param->s);

   if (invariant == CM_INVARIANT_SIMPLEETA)
      param->field = CM_FIELD_COMPLEX;
   else if (invariant == CM_INVARIANT_MULTIETA) {
      for (i = 0; param->p [i] != 0; i++);
      if (i % 2 == 0)
         param->field = CM_FIELD_REAL;
      else
         param->field = CM_FIELD_COMPLEX;
   }
   else
      param->field = CM_FIELD_REAL;

   return true;
}

/*****************************************************************************/

static bool doubleeta_compute_parameter (cm_param_ptr param, int_cl_t d)
   /* Compute p1 <= p2 prime following Cor. 3.1 of [EnSc04], that is,
      - 24 | (p1-1)(p2-1)
      - p1, p2 are not inert
      - if p1!=p2, then p1, p2 do not divide the conductor
      - if p1=p2=p, then either p splits or divides the conductor.
      The case p1=p2=2 of Cor. 3.1 is excluded by divisibility by 24.
      Minimise with respect to the height factor gained, which is
      12 psi (p1*p2) / (p1-1)(p2-1);
      then p1, p2 <= the smallest split prime which is 1 (mod 24). */

{
   int_cl_t cond2 = d / cm_classgroup_fundamental_discriminant (d);
      /* square of conductor */
   const unsigned long int maxprime = 997;
   unsigned long int primelist [169];
      /* list of suitable primes, terminated by 0; big enough to hold all
         primes <= maxprime */
   unsigned long int p1, p2, p1opt = 0, p2opt = 0;
   double quality, opt;
   bool ok;
   int i, j;

   /* determine all non-inert primes */
   i = 0;
   p1 = 2;
   ok = false;
   do {
      int kro = cm_nt_kronecker (d, (int_cl_t) p1);
      if (kro != -1) {
         primelist [i] = p1;
         i++;
      }
      if (kro == 1 && (p1 - 1) % 24 == 0)
         ok = true;
      else
         p1 = cm_nt_next_prime (p1);
   }
   while (p1 <= maxprime && !ok);
   primelist [i] = 0;

   /* search for the best tuple */
   opt = 0.0;
   for (j = 0, p2 = primelist [j]; p2 != 0; j++, p2 = primelist [j])
      for (i = 0, p1 = primelist [i]; i <= j; i++, p1 = primelist [i])
         if (   ((p1 - 1)*(p2 - 1)) % 24 == 0
             && (   (p1 != p2 && cond2 % p1 != 0 && cond2 % p2 != 0)
                 || (p1 == p2 && ((-d) % p1 != 0 || cond2 % p1 == 0)))) {
            quality = (p1 == p2 ? p1 : p1 + 1) * (p2 + 1) / (double) (p1 - 1)
               / (double) (p2 - 1);
            if (quality > opt) {
               p1opt = p1;
               p2opt = p2;
               opt = quality;
               ok = true;
            }
         }

   param->p [0] = p1opt;
   param->p [1] = p2opt;
   param->p [2] = 0;
   param->s = 1;
   param->e = 1;

   return ok;
}

/*****************************************************************************/
