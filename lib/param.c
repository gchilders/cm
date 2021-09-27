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

static void param_string (char *str, cm_param_ptr param);
static bool simpleeta_compute_parameter (cm_param_ptr param, int_cl_t d);
static void doubleeta_compute_parameter (cm_param_ptr param,
   cm_param_ptr paramsf, int_cl_t d, int maxdeg);
static void multieta_compute_parameter (cm_param_ptr param,
   cm_param_ptr paramsf2, cm_param_ptr paramsf4, int_cl_t d,
   int maxdeg);

/*****************************************************************************/

bool cm_param_init (cm_param_ptr param, int_cl_t d, char invariant,
   int maxdeg, int subfield, bool verbose)
   /* Test whether the discriminant is suited for the chosen invariant and
      in this case compute and store the parameters in param and return
      true; otherwise return false.
      If it is positive, then maxdeg is an upper bound on the degree of
      the modular polynomial in j, which is taken into account for certain
      infinite families of class invariants. If maxdeg is set to -1, then
      an internal bound is activated depending on the type of invariant.
      subfield is one of the following three constants:
      - CM_SUBFIELD_NEVER: Indicates that a generator of the class field
        is desired.
      - CM_SUBFIELD_PREFERRED: Indicates that whenever possible, a subfield
        of minimal degree of the class field should be computed, even if
        this leads to a worse height. This should be preferable for ECPP,
        where most of the time is spent in factoring class polynomials.
      - CM_SUBFIELD_OPTIMAL: Indicates that the smallest polynomial should
        be computed, taking the degree and the height of the class
        polynomial into account. */
{
   cm_param_t paramsf2, paramsf4;
   double size, sizesf2, sizesf4;

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
   param->r [0] = 0;
      /* Default choices; may be overwritten below. */
   memcpy (paramsf2, param, sizeof (__cm_param_struct));
   memcpy (paramsf4, param, sizeof (__cm_param_struct));

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
         doubleeta_compute_parameter (param, paramsf2, d, maxdeg);
         if (subfield == CM_SUBFIELD_NEVER) {
            if (param->p [0] == 0)
               return false;
         }
         else if (subfield == CM_SUBFIELD_PREFERRED) {
            if (paramsf2->p [0] != 0)
               param [0] = paramsf2 [0];
            else if (param->p [0] == 0)
               return false;
         }
         else /* CM_SUBFIELD_OPTIMAL */ {
            if (param->p [0] == 0)
               if (paramsf2->p [0] == 0)
                  return false;
               else
                  param [0] = paramsf2 [0];
            else if (paramsf2->p [0] != 0) {
               size = cm_class_height_factor (param);
               sizesf2 = 2 * cm_class_height_factor (paramsf2);
                  /* The computed height already takes into account that the
                     size of the coefficients is divided by about 2; the
                     additional factor 2 takes the degree halving into
                     account, so that the total size of the class polynomial
                     is measured. */
               if (sizesf2 > size)
                  param [0] = paramsf2 [0];
            }
         }
         break;
      case CM_INVARIANT_MULTIETA:
         multieta_compute_parameter (param, paramsf2, paramsf4, d, maxdeg);
         if (subfield == CM_SUBFIELD_NEVER) {
            if (param->p [0] == 0)
               return false;
         }
         else if (subfield == CM_SUBFIELD_PREFERRED) {
            if (paramsf4->p [0] != 0) {
               param [0] = paramsf4 [0];
               param->field = CM_FIELD_COMPLEX;
                  /* should not be touched, but identification of complex
                     conjugate values is not yet identified */
            }
            else if (paramsf2->p [0] != 0){
               param [0] = paramsf2 [0];
               param->field = CM_FIELD_COMPLEX; /* should not be touched */
            }
            else if (param->p [0] == 0)
               return false;
         }
         else /* CM_SUBFIELD_OPTIMAL */ {
            if (param->p [0] == 0 && paramsf2->p [0] == 0
                && paramsf4->p [0] == 0)
                return false;
            size = (param->p [0] == 0 ? 0.0 :
                    cm_class_height_factor (param));
            sizesf2 = (paramsf2->p [0] == 0 ? 0.0 :
                       2 * cm_class_height_factor (paramsf2));
            sizesf4 = (paramsf4->p [0] == 0 ? 0.0 :
                       4 * cm_class_height_factor (paramsf4));
               /* The height factor already takes the reduced height due
                  to subfields into account, the additional factors 2 and 4
                  correspond to the gain in the number of factors. */
            if (sizesf4 >= sizesf2 && sizesf4 >= size) {
               param [0] = paramsf4 [0];
               param->field = CM_FIELD_COMPLEX; /* should not be touched */
            }
            else if (sizesf2 >= size) {
               param [0] = paramsf2 [0];
               param->field = CM_FIELD_COMPLEX; /* should not be touched */
            }
         }
         break;
      default: /* should not occur */
         printf ("class_compute_parameter called for "
                  "unknown class invariant '%c'\n", invariant);
         exit (1);
   }

   /* Create parameter string. */
   param_string (param->str, param);

   return true;
}

/*****************************************************************************/

static void param_string (char *str, cm_param_ptr param)
   /* Write the numerical parameters from param into str, which needs to be
      allocated with sufficient space. */
{
   int i;

   if (param->p [0] != 0) {
      i = 0;
      do {
         str += sprintf (str, "%i_", param->p [i]);
         i++;
      } while (param->p [i] != 0);
      str += sprintf (str, "%i_%i", param->e, param->s);
      for (i = 0; param->r [i] != 0; i++)
         str += sprintf (str, "_%i", param->r [i]);
   }
   else
      str [0] = '\0';
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

static void doubleeta_compute_parameter (cm_param_ptr param,
   cm_param_ptr paramsf, int_cl_t d, int maxdeg)
   /* Compute and return in param a parameter combination for a double eta
      quotient that, as far as we know, does not lead to a subfield; and in
      paramsf a parameter combination that leads to a subfield of index 2.
      maxdeg has the same meaning as in cm_param_init; if set to -1, it is
      internally replaced by 2, in which case the modular curve
      X_0^+ (p1*p2) has genus 0; otherwise said, both roots of the modular
      polynomial in j lead to curves of the correct cardinality.
      If a suitable param or paramsf does not exist, this is indicated
      by setting their p [0] value to 0.
      It is up to the calling function to decide which of param or paramsf
      to use.
      Compute p1 <= p2 prime and s following Cor. 3.1 of [EnSc04], that is,
      - s = 24 / gcd (24, (p1-1)(p2-1))
      - p1, p2 are not inert
      - if p1!=p2, then p1, p2 do not divide the conductor
      - if p1=p2=p!=2, then either p splits or divides the conductor
      - if p1=p2=2, then either 2 splits, or 2 divides the conductor
        and d != 4 (mod 32).
      Additionally consider the lower powers e given in Theorem 1 of
      [EnSc13] for p1 != p2, and none of them inert or dividing the
      conductor.
      For the subfield case, use Theorem 5 of [EnSc13]. A subfield of
      index 2 is obtained if additionally to the previous condition,
      p1!=p2 are ramified and d is neither -p1*p2 nor -4*p1*p2, and
      one of the following condition holds:
      - e is even;
      - p1, p2 != 2 and at least one of them is 1 mod 4;
      - p1 = 2 and p2 = +- 1 mod 8.
      If none of these three conditions holds for the optimal e and
      s is even (which implies that p1 and p2 are 3 mod 4), then we
      may also use 2*e.
      Minimise with respect to the height factor gained. */
{
   int_cl_t cond2 = d / cm_classgroup_fundamental_discriminant (d);
      /* square of conductor */
   const long int maxprime = 997;
   long int primelist [168];
      /* list of suitable primes; big enough to hold all
         primes <= maxprime */
   int length; /* effective length of primelist */
   int p, p1, p2, s, e;
   cm_param_t par;
   double opt, optsf;
   double hf;
   int i, j;

   if (maxdeg == -1)
      maxdeg = 2;

   /* Determine all non-inert primes up to maxprime. */
   length = 0;
   for (p = 2; p < maxprime; p = cm_nt_next_prime (p))
      if (cm_nt_kronecker (d, (int_cl_t) p) != -1) {
         primelist [length] = p;
         length++;
      }

   /* Search for the best tuple. */
   opt = 0.0;
   optsf = 0.0;
   param->p [0] = 0;
   paramsf->p [0] = 0;
   par [0] = param [0]; /* copy d and invariant fields */
   par->p [2] = 0;
   par->r [2] = 0;
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

            /* Choose e according to [EnSc13], Theorem 1. */
            if (p1 != p2 && p1 != 2) {
               if (p1 == 3 || d % 3 == 0)
                  e = 3 / cm_nt_gcd (3, (p1 - 1) * (p2 - 1));
               else
                  e = 1;
               if (d % 2 == 0)
                  e *= 8 / cm_nt_gcd (8, (p1 - 1) * (p2 - 1));
            }
            else
               e = s;

            par->p [0] = p1;
            par->p [1] = p2;
            par->s = s;
            par->e = e;
            hf = cm_class_height_factor (par);
            if (hf > opt || hf > optsf) {
               if (p1 == p2 || d % p1 != 0 || d % p2 != 0
                   || d == -p1*p2 || d == -4*p1*p2) {
                  /* We are not in the subfield case regardless of e. */
                  if (hf > opt) {
                     param [0] = par [0];
                     opt = hf;
                  }
               }
               else if (e % 2 == 0
                        || (p1 != 2 && (p1 % 4 == 1 || p2 % 4 == 1))
                        || (p1 == 2 && (p2 % 8 == 1 || p2 % 8 == 7))) {
                  /* We are in the subfield case. */
                  if (hf > optsf) {
                     paramsf [0] = par [0];
                     paramsf->r [0] = p1;
                     paramsf->r [1] = p2;
                     optsf = hf;
                  }
               }
               else if (s % 2 == 0) {
                  /* We are in the subfield case for 2*e. */
                  hf /= 2;
                  par->e *= 2;
                  if (hf > optsf) {
                     paramsf [0] = par [0];
                     paramsf->r [0] = p1;
                     paramsf->r [1] = p2;
                     optsf = hf;
                  }
               }
            }
         }
      }
   }
}

/*****************************************************************************/

static void multieta_compute_parameter (cm_param_ptr param,
   cm_param_ptr paramsf2, cm_param_ptr paramsf4, int_cl_t d, int maxdeg)
   /* Choose parameter combinations for triple eta quotients that (as far
      as we know) lead to a full class field, a subfield of index 2 in the
      class field or a subfield of index 4. If no suitable choice exists,
      then the p [0] value of the corresponding parameter is set to 0.
      We consider only the multiple eta quotients for which the modular
      polynomial has a degree of at most 8 in j. These are:
      p1,p2,p3  s  hf   deg
      2,3,5     3  18    4
      2,3,7     2  24    4
      2,3,13    1  42    4
      2,5,7     1  36    4
      2,5,13    1  31.5  8
      3,5,7     1  24    8
      Of interest could also be the smallest case with four factors, but
      even in the totally ramified case its degree in j is not optimal:
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
   int_cl_t dfund = cm_classgroup_fundamental_discriminant (d);
   int_cl_t cond2 = d / dfund;
   int prime [] = {2, 3, 5, 7, 13};
   const int prime_l = sizeof (prime) / sizeof (int);
   cm_param_t cand [6];
   int deg [6];
   const int cand_l = sizeof (cand) / sizeof (cm_param_t);
   bool ok [14], ram [14];
   int i, j;

   if (maxdeg == -1)
      maxdeg = 4;

   for (i = 0; i < prime_l; i++) {
      ok [prime [i]] = cm_nt_kronecker (d, (int_cl_t) prime [i]) != -1
                       && cond2 % prime [i] != 0;
      ram [prime [i]] = (dfund % prime [i] == 0);
   }

   /* Initialise parameter combinations with decreasing height factors
      and, for the same height factor, with increasing degree in j. */
   for (i = 0; i < cand_l; i++) {
      cand [i][0] = param [0];
         /* Copies in particular the discriminant. */
      cand [i]->p [0] = 2;
      cand [i]->p [1] = 3;
      cand [i]->p [2] = 13;
      cand [i]->p [3] = 0;
      cand [i]->s = 1;
      cand [i]->e = 1;
      cand [i]->r [0] = 0;
      deg [i] = 4;
   }
   cand [1]->p [1] = 5;
   cand [1]->p [2] = 7;
   cand [2]->p [1] = 5;
   deg [2] = 8;
   cand [3]->p [2] = 7;
   cand [3]->s = 2;
   cand [3]->e = 2;
   cand [4]->p [0] = 3;
   cand [4]->p [1] = 5;
   cand [4]->p [2] = 7;
   deg [4] = 8;
   cand [5]->p [2] = 5;
   cand [5]->s = 3;
   cand [5]->e = 3;

   param->p [0] = 0;
   paramsf2->p [0] = 0;
   paramsf4->p [0] = 0;
   for (i = 0; i < cand_l; i++) {
      if (   (maxdeg == 0 || deg [i] <= maxdeg)
          && ok [cand [i]->p [0]]
          && ok [cand [i]->p [1]]
          && ok [cand [i]->p [2]]) {
         /* Parameter combination yields a class invariant, check
            ramification and determine subfield. */
         if (ram [cand [i]->p [0]] && ram [cand [i]->p [1]]
            && ram [cand [i]->p [2]]) {
            if (paramsf4->p [0] == 0) {
               paramsf4 [0] = cand [i][0];
               for (j = 0; j < 3; j++)
                  paramsf4->r [j] = cand [i]->p [j];
            }
         }
         else if (ram [cand [i]->p [0]] && ram [cand [i]->p [1]]) {
            if (paramsf2->p [0] == 0) {
               paramsf2 [0] = cand [i][0];
               paramsf2->r [0] = cand [i]->p [0];
               paramsf2->r [1] = cand [i]->p [1];
            }
         }
         else if (ram [cand [i]->p [0]] && ram [cand [i]->p [2]]) {
            if (paramsf2->p [0] == 0) {
               paramsf2 [0] = cand [i][0];
               paramsf2->r [0] = cand [i]->p [0];
               paramsf2->r [1] = cand [i]->p [2];
            }
         }
         else if (ram [cand [i]->p [1]] && ram [cand [i]->p [2]]) {
            if (paramsf2->p [0] == 0) {
               paramsf2 [0] = cand [i][0];
               paramsf2->r [0] = cand [i]->p [1];
               paramsf2->r [1] = cand [i]->p [2];
            }
         }
         else {
            if (param->p [0] == 0)
               param [0] = cand [i][0];
         }
      }
   }
   paramsf2->r [2] = 0;
   paramsf4->r [3] = 0;

   /* The polynomial is real by [EnSc13] if the number of primes is
      even (Corollary 4), or any of the primes is ramified (Corollary 8). */
   param->field = CM_FIELD_COMPLEX;
   if (param->p [0] != 0)
      for (j = 0; j < 3; j++)
         if (ram [param->p [j]])
            param->field = CM_FIELD_REAL;
   paramsf2->field = CM_FIELD_REAL;
   paramsf4->field = CM_FIELD_REAL;
}

/*****************************************************************************/
