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

static bool simpleeta_compute_parameter (cm_param_ptr param, int_cl_t d);
static bool doubleeta_compute_parameter (cm_param_ptr param, int_cl_t d,
   int maxdeg);
static bool multieta_compute_parameter (cm_param_ptr param, int_cl_t d,
   int maxdeg);

/*****************************************************************************/

bool cm_param_init (cm_param_ptr param, int_cl_t d, char invariant,
   int maxdeg, bool verbose)
   /* Test whether the discriminant is suited for the chosen invariant and
      in this case compute and store the parameters in param and return
      true; otherwise return false.
      If it is positive, then maxdeg is an upper bound on the degree of
      the modular polynomial in j, which is taken into account for certain
      infinite families of class invariants. If maxdeg is set to -1, then
      an internal bound is activated depending on the type of invariant. */

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
   param->field = CM_FIELD_REAL;
      /* Default choice; may be overwritten below. */

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
      case CM_INVARIANT_SIMPLEETA:
         if (!simpleeta_compute_parameter (param, d))
            return false;
         break;
      case CM_INVARIANT_DOUBLEETA:
         if (!doubleeta_compute_parameter (param, d, maxdeg))
            return false;
         break;
      case CM_INVARIANT_MULTIETA:
         if (!multieta_compute_parameter (param, d, maxdeg))
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

   return true;
}

/*****************************************************************************/

static bool simpleeta_compute_parameter (cm_param_ptr param, int_cl_t d)
   /* If any exist, compute n and e following [EnMo14] such that w_n^e is a
      class invariant for d and j can be obtained from w_n without the use
      of modular polynomials (otherwise said, n is one of 3, 5, 7, 13,
      4, 9 or 25). If several exist, return in param the one with the best
      height factor. */
{
   int k3, k5, k7, k13;
   uint_cl_t dmod12, dmod36, dmod108, dmod8, dmod32, dmod64, dmod128;

   k3 = cm_nt_kronecker (d, (int_cl_t) 3);
   k5 = cm_nt_kronecker (d, (int_cl_t) 5);
   k7 = cm_nt_kronecker (d, (int_cl_t) 7);
   k13 = cm_nt_kronecker (d, (int_cl_t) 13);
   dmod12 = cm_classgroup_mod (d, (uint_cl_t) 12);
   dmod36 = cm_classgroup_mod (d, (uint_cl_t) 36);
   dmod108 = cm_classgroup_mod (d, (uint_cl_t) 108);
   dmod8 = cm_classgroup_mod (d, (uint_cl_t) 8);
   dmod32 = cm_classgroup_mod (d, (uint_cl_t) 32);
   dmod64 = cm_classgroup_mod (d, (uint_cl_t) 64);
   dmod128 = cm_classgroup_mod (d, (uint_cl_t) 128);

   /* According to [EnMo04], we have s=24/(n-1), and the height factor is
      l+1 for n=l prime and l for n=l^2. Moreover, w_l^s can be used
      whenever kronecker(d,l) != -1, and w^{l^2}^s can be used whenever
      kronecker(d,l) == 1 and for some non-fundamental discriminants with
      kronecker(d,l) == 0; so full powers of squares other than 4 are not
      of interest, and in particular 25 can be dropped.
      Test all possible powers in the order of their height factors. */
   if (k3 != -1 &&
       (dmod12 == 1 || dmod36 == 33)) {
      /* w_3^2, factor 24 */
      param->p [0] = 3;
      param->e = 2;
   }
   else if (k5 != -1 && d % 3 != 0) {
      /* w_5^2, factor 18 */
      param->p [0] = 5;
      param->e = 2;
   }
   else if (k7 != -1 && d % 2 != 0) {
      /* w_7^2, factor 16 */
      param->p [0] = 7;
      param->e = 2;
   }
   else if (dmod8 == 1 || dmod64 == 48 || dmod128 == 0) {
      /* w_4, factor 16 */
      param->p [0] = 4;
      param->e = 1;
   }
   else if (k13 != -1) {
      /* w_13^2, factor 14 */
      param->p [0] = 13;
      param->e = 2;
   }
   else if (k3 != -1 &&
       (dmod12 == 4 || dmod36 == 24)) {
      /* w_3^4, factor 12 */
      param->p [0] = 3;
      param->e = 4;
   }
   else if (   dmod108 == 0 || dmod108 == 45
            || dmod108 == 72 || dmod108 == 81) {
      /* w_9, factor 9 */
      param->p [0] = 9;
      param->e = 1;
   }
   else if (k7 != -1) {
      /* d even, w_7^4, factor 8 */
      param->p [0] = 7;
      param->e = 4;
   }
   else if (k3 != -1 &&
       (dmod36 == 9 || dmod36 == 21)) {
      /* w_3^6, factor 8 */
      param->p [0] = 3;
      param->e = 6;
   }
   else if (dmod32 == 20 || dmod128 == 64) {
      /* w_4^2, factor 8 */
      param->p [0] = 4;
      param->e = 2;
   }
   else if (k5 != -1) {
      /* 3|d, w_5^6, factor 6 */
      param->p [0] = 5;
      param->e = 6;
   }
   else if (k3 != -1) {
      /* w_3^12, factor 4 */
      param->p [0] = 3;
      param->e = 12;
   }
   else if (dmod64 == 16 || dmod64 == 32) {
      /* w_4^4, factor 4 */
      param->p [0] = 4;
      param->e = 4;
   }
   else if (dmod8 == 1 || dmod32 == 4) {
      /* w4^8, factor 2 */
      param->p [0] = 4;
      param->e = 8;
   }
   else
      return false;

   param->p [1] = 0;
   param->s = 24 / (param->p [0] - 1);

   /* The field is complex by default, but real in some cases worked out
      in [EnMo14], Theorems 4.4 and 6.1. The condition 16|d for p [0] == 4
      in Theorem 6.1 is only necessary, but not sufficient. Looking
      more closely at the conditions for N-systems shows that for
      p [0] == 4, d = 16 (mod 64) is the only case in which our choice of
      class invariants above provenly yields a real class polynomial.
      (Whenever 16|d and s==e, one also obtains a real polynomial; but then
      one can always choose a lower power, which is usually not real.) */
   if (param->p [0] == 4)
      if (dmod64 == 16)
         param->field = CM_FIELD_REAL;
      else
         param->field = CM_FIELD_COMPLEX;
   else
      if (d % param->p [0] == 0
          && (param->e == param->s
              || (param->p [0] == 3 && param->e == 4
                  && cm_classgroup_mod (d, 9) == 6)
              || (param->p [0] == 5 && d % 3 != 0)
              || (param->p [0] == 13 && cm_classgroup_mod (d, 27) == 18)))
         param->field = CM_FIELD_REAL;
      else
         param->field = CM_FIELD_COMPLEX;

   return true;
}

/*****************************************************************************/

static bool doubleeta_compute_parameter (cm_param_ptr param, int_cl_t d,
   int maxdeg)
   /* Compute p1 <= p2 prime and s following Cor. 3.1 of [EnSc04], that is,
      - s = 24 / gcd (24, (p1-1)(p2-1))
      - p1, p2 are not inert
      - if p1!=p2, then p1, p2 do not divide the conductor
      - if p1=p2=p!=2, then either p splits or divides the conductor
      - if p1=p2=2, then either 2 splits, or 2 divides the conductor
        and d != 4 (mod 32).
      Additionally consider the lower powers e given in Theorem 1 of
      [EnSc13] for p1 != p2, and none of them inert or dividing the
      conductor.
      Minimise with respect to the height factor gained; then the
      optimal primes are bounded above by the smallest split prime
      that is 1 (mod 24).
      maxdeg has the same meaning as in cm_param_init; if set to -1, it is
      internally replaced by 2, in which case the modular curve
      X_0^+ (p1*p2) has genus 0; otherwise said, both roots of the modular
      polynomial in j lead to curves of the correct cardinality. */

{
   int_cl_t cond2 = d / cm_classgroup_fundamental_discriminant (d);
      /* square of conductor */
   const long int maxprime = 997;
   long int primelist [168];
      /* list of suitable primes; big enough to hold all
         primes <= maxprime */
   int length; /* effective length of primelist */
   long int p, p1, p2, s;
   cm_param_t par;
   double hf, opt;
   bool ok;
   int i, j;

   if (maxdeg == -1)
      maxdeg = 2;

   /* Determine all non-inert primes up to maxprime or a split prime that
      is 1 modulo 24, whichever comes first. */
   length = 0;
   p = 2;
   ok = false;
   do {
      int kro = cm_nt_kronecker (d, (int_cl_t) p);
      if (kro != -1) {
         primelist [length] = p;
         length++;
      }
      if (kro == 1 && (p - 1) % 24 == 0)
         ok = true;
      else
         p = cm_nt_next_prime (p);
   }
   while (p <= maxprime && !ok);

   /* Search for the best tuple. */
   opt = 0.0;
   par [0] = param [0];
   par->p [2] = 0;
   ok = false;
   for (j = 0; j < length; j++) {
      p2 = primelist [j];
      for (i = 0; i <= j; i++) {
         p1 = primelist [i];
         s = 24 / cm_nt_gcd (24, (p1 - 1) * (p2 - 1));
         if (((p1 != p2 && cond2 % p1 != 0 && cond2 % p2 != 0)
               || (p1 == p2 && p1 != 2
                   && (d % p1 != 0 || cond2 % p1 == 0))
               || (p1 == 2 && p2 == 2
                   && (d % 2 != 0
                       || (cond2 % 2 == 0
                           && cm_classgroup_mod (d, 32) != 4))))
             && (maxdeg == 0
                 || s * (p1 - 1) * (p2 - 1) / 12 <= maxdeg)) {
            par->p [0] = p1;
            par->p [1] = p2;
            par->s = s;
            /* Choose e according to [EnSc13], Theorem 1. */
            if (p1 != p2 && p1 != 2) {
               if (p1 == 3 || d % 3 == 0)
                  par->e = 3 / cm_nt_gcd (3, (p1 - 1) * (p2 - 1));
               else
                  par->e = 1;
               if (d % 2 == 0)
                  par->e *= 8 / cm_nt_gcd (8, (p1 - 1) * (p2 - 1));
            }
            else
               par->e = s;
            hf = cm_class_height_factor (par);
            if (hf > opt) {
               param [0] = par [0];
               opt = hf;
               ok = true;
            }
         }
      }
   }

   return ok;
}

/*****************************************************************************/

static bool multieta_compute_parameter (cm_param_ptr param, int_cl_t d,
   int maxdeg)
   /* We only consider the multiple eta quotients for which the modular
      polynomial has a degree of at most 8 in j. These are:
      p1,p2,p3  s  hf   deg
      2,3,5     3  18    4
      2,3,7     2  24    4
      2,3,13    1  42    4
      2,5,7     1  36    4
      2,5,13    1  31.5  8
      3,5,7     1  24    8
      Of interest could also be the smallest case with four factors,
      but even in the ramified case its degree in j is not optimal:
      2,3,5,7   1  36   16
      If one wanted to go up to a degree 16, the following quotients
      with three primes would have to be added:
      2,3,19    2  20   12
      2,3,37    1  38   12
      2,5,19    1  30   12
      2,7,13    1  28   12
      2,3,17    3  13.5 16
      2,7,17    1  27   16
      3,5,13    1  21   16 */
{
   int_cl_t cond2 = d / cm_classgroup_fundamental_discriminant (d);
      /* square of conductor */
   bool ok2, ok3, ok5, ok7, ok13, ramified;
   int k;

   if (maxdeg == -1 || (maxdeg >= 4 && maxdeg < 8))
      maxdeg == 4;
   else if (maxdeg > 0 && maxdeg < 4)
      return false;
   else
      maxdeg = 8;

   ok2 = cm_nt_kronecker (d, (int_cl_t) 2) != -1 && cond2 % 2 != 0;
   ok3 = cm_nt_kronecker (d, (int_cl_t) 3) != -1 && cond2 % 3 != 0;
   ok5 = cm_nt_kronecker (d, (int_cl_t) 5) != -1 && cond2 % 5 != 0;
   ok7 = cm_nt_kronecker (d, (int_cl_t) 7) != -1 && cond2 % 7 != 0;
   ok13 = cm_nt_kronecker (d, (int_cl_t) 13) != -1 && cond2 % 13 != 0;

   param->p [3] = 0;
   param->s = 1; /* default, can be overwritten */
   if (ok2) {
      param->p [0] = 2;
      if (ok3 && ok13) {
         /* height factor 42 */
         param->p [1] = 3;
         param->p [2] = 13;
      }
      else if (ok5 && ok7) {
         /* height factor 36 */
         param->p [1] = 5;
         param->p [2] = 7;
      }
      else if (maxdeg >= 8 && ok5 && ok13) {
         /* height factor 31.5 */
         param->p [1] = 5;
         param->p [2] = 13;
      }
      else if (ok3 && ok7) {
         /* height factor 24 */
         param->p [1] = 3;
         param->p [2] = 7;
         param->s = 2;
      }
      else if (ok3 && ok5) {
         /* height factor 18 */
         param->p [1] = 3;
         param->p [2] = 5;
         param->s = 3;
      }
      else
         return false;
   }
   else if (maxdeg >= 8 && ok3 && ok5 && ok7) {
      param->p [0] = 3;
      param->p [1] = 5;
      param->p [2] = 7;
   }
   else
      return false;
   param->e = param->s;

   /* The polynomial is real by [EnSc13] if the number of primes is
      even (Corollary 4), or any of the primes is ramified (Corollary 8).
      Currently we use exactly three prime factors, but there is no harm
      in writing more complete code already now. */
   ramified = false;
   for (k = 0; param->p [k] != 0; k++)
      if (!ramified && d % param->p [k] == 0)
         ramified = true;
   if (k % 2 == 0 || ramified)
      param->field = CM_FIELD_REAL;
   else
      param->field = CM_FIELD_COMPLEX;

   return true;
}

/*****************************************************************************/
