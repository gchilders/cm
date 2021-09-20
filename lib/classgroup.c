/*

classgroup.c - computations with class groups and quadratic forms

Copyright (C) 2009, 2010, 2018, 2021 Andreas Enge

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

static int_cl_t classgroup_gcdext (int_cl_t *u, int_cl_t *v, int_cl_t a,
   int_cl_t b);

/*****************************************************************************/
/*                                                                           */
/* functions transforming between int_cl_t and mpz_t                         */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_mpz_set_icl (mpz_t rop, int_cl_t op)
   /* relies on the fact that int_cl_t has 64 bits */

{
   int_cl_t mask = (((int_cl_t) 1) << 32) - 1;
   int_cl_t abso = (op > 0 ? op : -op);

   /* copy 32 high bits */
   mpz_set_ui (rop, (unsigned long int) (abso >> 32));
   /* add 32 low bits */
   mpz_mul_2exp (rop, rop, 32);
   mpz_add_ui (rop, rop, (unsigned long int) (abso & mask));

   if (op < 0)
      mpz_neg (rop, rop);
}

/*****************************************************************************/

int_cl_t cm_classgroup_mpz_get_icl (mpz_t op)

{
   int_cl_t rop;
   int sign = (mpz_cmp_ui (op, 0) < 0 ? -1 : 1);
   mpz_t tmp;

   mpz_init (tmp);

   /* switch to absolute value */
   if (sign < 0)
      mpz_neg (op, op);

   /* extract 32 high bits */
   mpz_tdiv_q_2exp (tmp, op, 32);
   rop = ((int_cl_t) (mpz_get_ui (tmp))) << 32;
   /* extract 32 low bits */
   mpz_tdiv_r_2exp (tmp, op, 32);
   rop += mpz_get_ui (tmp);

   if (sign < 0) {
      rop = -rop;
      mpz_neg (op, op);
   }

   mpz_clear (tmp);

   return rop;
}


/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_init (cm_classgroup_t *cl, int_cl_t disc, bool verbose)

{
   int h;
   int length;
   /* Nobody will ever need a class group with more than 64 components.
      More precisely, this is the limit imposed by the class number encoded
      in 64 bits. */
   int_cl_t ord [64];
   cm_form_t gen [64];
   cm_form_t Ppow;
   int i, j, k;

   if (disc >= 0) {
      printf ("\n*** The discriminant must be negative.\n");
      exit (1);
   }
   else if (disc % 4 != 0 && (disc - 1) % 4 != 0) {
      printf ("\n*** %"PRIicl" is not a quadratic discriminant.\n", disc);
      exit (1);
   }
   else
      cl->d = disc;

   length = cm_classgroup_normalseries (disc, NULL, ord, gen);
   cl->h = 1;
   for (i = 0; i < length; i++)
      cl->h *= ord [i];
   cl->form = (cm_form_t *) malloc (cl->h * sizeof (cm_form_t));
   cl->conj = (int *) malloc (cl->h * sizeof (int));
   if (verbose)
      printf ("Class number: h = %d\n", cl->h);

   cl->levels = length;
   cl->deg = (int *) malloc (length * sizeof (int));
   for (i = 0; i < length; i++)
      cl->deg [i] = ord [length - 1 - i];
   if (verbose)
      for (i = 0; i < length; i++)
         printf ("%"PRIicl" [%"PRIicl", %"PRIicl"]\n",
                 ord [i], gen [i].a, gen [i].b);

   /* Roll out the class group from its generators. */
   cl->form [0].a = 1;
   cl->form [0].b = (disc % 2 == 0 ? 0 : 1);
   h = 1;
   for (i = 0; i < length; i++) {
      Ppow = gen [i];
      for (j = 1; j < ord [i]; j++) {
         /* Ppow contains gen [i]^j; multiply the first h forms by this. */
         for (k = 0; k < h; k++)
            cm_classgroup_compose (&(cl->form [h*j+k]), cl->form [k], Ppow,
                                   disc);
         cm_classgroup_compose (&Ppow, Ppow, gen [i], disc);
      }
      h *= ord [i];
   }

   /* Pair up inverse forms. */
   for (i = 0; i < cl->h; i++) {
      for (j = i;
           j < cl->h && cl->form [j].b != -cl->form [i].b;
           j++);
      if (j == cl->h || cl->form [i].b == 0)
         cl->conj [i] = i;
      else {
         cl->conj [i] = j;
         cl->conj [j] = i;
      }
   }
}

/*****************************************************************************/

void cm_classgroup_clear (cm_classgroup_t *cl)

{
   free (cl->form);
   free (cl->conj);
   free (cl->deg);
}

/*****************************************************************************/
/*                                                                           */
/* elementary number theory with int_cl_t                                    */
/*                                                                           */
/*****************************************************************************/

uint_cl_t cm_classgroup_mod (int_cl_t a, uint_cl_t p)
      /* returns a representative of a % p in [0, p-1[ */

{
   int_cl_t res = a % (int_cl_t) p;
   return (res < 0 ? res + (int_cl_t) p : res);
}

/*****************************************************************************/

int_cl_t cm_classgroup_gcd (int_cl_t a, int_cl_t b)
   /* returns the positive gcd of a and b */

{
   int_cl_t r;

   if (a == 0)
      return b;
   else if (b == 0)
      return a;
   else
   {
      if (a < 0)
         a = -a;
      if (b < 0)
         b = -b;
      r = a % b;
      while (r > 0)
      {
         a = b;
         b = r;
         r = a % b;
      }
      return b;
   }
}

/*****************************************************************************/

static int_cl_t classgroup_gcdext (int_cl_t *u, int_cl_t *v, int_cl_t a,
   int_cl_t b)
   /* returns the positive gcd d of a and b; if u and v is not NULL,         */
   /* modifies them such that u a + v b = d; it is also possible to have     */
   /* only one of u and v non NULL.                                          */

{
   int_cl_t r0, r1, r2, u0, u1, u2, v0, v1, v2, q;
   int sgn_a, sgn_b;

   if (a < 0) {
      sgn_a = -1;
      r0 = -a;
   }
   else {
      sgn_a = 1;
      r0 = a;
   }
   if (b < 0) {
      sgn_b = -1;
      r1 = -b;
   }
   else {
      sgn_b = 1;
      r1 = b;
   }
   /* computing u and v might be faster than permanent tests for NULL */
   u0 = 1;
   u1 = 0;
   v0 = 0;
   v1 = 1;

   while (r1 != 0) {
      q = r0 / r1;
      r2 = r0 % r1;
      r0 = r1;
      r1 = r2;
      u2 = u0 - q * u1;
      u0 = u1;
      u1 = u2;
      v2 = v0 - q * v1;
      v0 = v1;
      v1 = v2;
   }

   if (u != NULL)
      *u = sgn_a * u0;
   if (v != NULL)
      *v = sgn_b * v0;

   return r0;
}


/*****************************************************************************/
/*                                                                           */
/* functions computing fundamental discriminants, conductors and class       */
/* numbers                                                                   */
/*                                                                           */
/*****************************************************************************/

int_cl_t cm_classgroup_fundamental_primes (int_cl_t *primes, int_cl_t d)
   /* Given a discriminant d, compute and return its fundamental
      discriminant. At the same time, put the "signed primes" dividing the
      fundamental discriminant as a zero-terminated array into primes,
      which needs to be initialised with at least 17 entries. The "signed
      primes" are -4, 8, -8, p with p=1 mod 4 and -p with p=3 mod 4. */
{
   uint_cl_t p [17];
   unsigned int e [17];
   int i, j;

   /* Make the discriminant "2-fundamental", assuming it was a discriminant
      in the first place. */
   while (d % 4 == 0)
      d /= 4;
   if ((d - 1) % 4 != 0)
      d *= 4;

   cm_nt_factor (-d, p, e);

   i = 0; /* counter for p */
   j = 0; /* counter for primes */
   /* Handle the 2-part first; the difficulty is to decide between
      +8 and -8. */
   if (p [i] == 2) {
      if (e [i] == 2)
         primes [j] = -4;
      else if ((d/(-8) - 1) % 4 == 0)
         primes [j] = -8;
      else
         primes [j] = 8;
      i++;
      j++;
   }
   /* Now loop through the odd primes. */
   while (p [i] != 0) {
      if (e [i] % 2 == 1) {
         if ((p [i] - 1) % 4 == 0)
            primes [j] = p [i];
         else
            primes [j] = -p [i];
         j++;
      }
      i++;
   }

   primes [j] = 0;

   d = primes [0];
   for (i = 1; primes [i] != 0; i++)
      d *= primes [i];

   return d;
}

/*****************************************************************************/

int_cl_t cm_classgroup_fundamental_discriminant (int_cl_t d)
   /* Return the fundamental discriminant associated with the
      discriminant d. */
{
   int_cl_t p [17];

   return cm_classgroup_fundamental_primes (p, d);
}


/*****************************************************************************/
/*                                                                           */
/* functions for computing in the class group                                */
/*                                                                           */
/*****************************************************************************/

int_cl_t cm_classgroup_compute_c (int_cl_t a, int_cl_t b, int_cl_t d)
   /* computes c = (b^2 - d) / (4*a), potentially switching to multi-        */
   /* precision to avoid intermediate overflow                               */

{
   int_cl_t c;

   if ((b > 0 ? b : -b)
         < ((int_cl_t) 1) << (4 * sizeof (int_cl_t) - 2))
      c = (b * b - d) / a / 4;
   else {
      /* should happen only rarely */
      mpz_t az, dz, cz;
      mpz_init (az);
      mpz_init (cz);
      mpz_init (dz);
      cm_classgroup_mpz_set_icl (dz, d);
      cm_classgroup_mpz_set_icl (az, a);
      cm_classgroup_mpz_set_icl (cz, b);
      mpz_mul (cz, cz, cz);
      mpz_sub (cz, cz, dz);
      mpz_div (cz, cz, az);
      mpz_div_2exp (cz, cz, 2);
      c = cm_classgroup_mpz_get_icl (cz);
      mpz_clear (az);
      mpz_clear (cz);
      mpz_clear (dz);
   }

   return c;
}

/*****************************************************************************/

void cm_classgroup_reduce (cm_form_t *Q, int_cl_t d)
   /* reduces the quadratic form Q without checking if it belongs indeed to */
   /* the discriminant d.                                                   */

{
   int_cl_t c, a_minus_b, two_a, offset;
   bool reduced;

   reduced = false;
   while (!reduced){
      /* first step: obtain |b| <= a */
      a_minus_b = Q->a - Q->b;
      two_a = 2 * Q->a;
      if (a_minus_b < 0) {
         a_minus_b++;
         /* a trick to obtain the correct rounding */
         offset = a_minus_b / two_a;
         /* Since a_minus_b is negative, a negative remainder is computed, */
         /* so the quotient is effectively rounded to the nearest integer  */
         offset--;
         /* offset is (a-b) / (2a) floored */
         offset *= two_a;
         Q->b += offset;
      }
      else if (a_minus_b >= two_a) {
         offset = a_minus_b / two_a;
         offset *= two_a;
         Q->b += offset;
      }

      c = cm_classgroup_compute_c (Q->a, Q->b, d);

      /* if not reduced, invert */
      if (Q->a < c || (Q->a == c && Q->b >= 0))
            reduced = true;
      else {
         Q->a = c;
         Q->b = -Q->b;
      }
   }
}

/*****************************************************************************/

void cm_classgroup_compose (cm_form_t *Q, cm_form_t Q1, cm_form_t Q2,
   int_cl_t d)
   /* computes the reduced form Q corresponding to the composition of Q1 and */
   /* Q2.                                                                     */

{
   int_cl_t s, t, v1, v, w, a2t;

   t = classgroup_gcdext (&v1, NULL, Q2.a, Q1.a);

   if (t == 1) {
      Q->a = Q1.a * Q2.a;
      Q->b = Q2.b + Q2.a * v1 * (Q1.b - Q2.b);
   }
   else {
      s = (Q1.b + Q2.b) / 2;
      t = classgroup_gcdext (&w, &v, s, t);
      v *= v1;
      a2t = Q2.a / t;
      Q->a = (Q1.a / t) * a2t;
      Q->b = ((s - Q2.b) * v - w * (Q2.b * Q2.b - d) / (4 * Q2.a))
             % (2 * Q->a); /* intermediate reduction to avoid overflow */
      Q->b = (Q2.b + 2 * Q->b * a2t) % (2 * Q->a);
   }
   cm_classgroup_reduce (Q, d);
}

/*****************************************************************************/

void cm_classgroup_pow (cm_form_t *Q, cm_form_t P, uint_cl_t n, int_cl_t d)
   /* Compute the reduced form Q corresponding to the n-th power of P. */

{
   cm_form_t R;

   if (n == 0) {
      Q->a = 1;
      Q->b = cm_classgroup_mod (d, 2);
   }
   else if (n % 2 == 0) {
      cm_classgroup_pow (&R, P, n/2, d);
      cm_classgroup_compose (Q, R, R, d);
         /* For more efficiency, which is not an issue here, we could use
            a dedicated squaring function. */
   }
   else /* n odd */ {
      cm_classgroup_pow (&R, P, n-1, d);
      cm_classgroup_compose (Q, R, P, d);
   }
}

/*****************************************************************************/
/*                                                                           */
/* function for computing a normal series of the class group                 */
/*                                                                           */
/*****************************************************************************/

int cm_classgroup_normalseries (int_cl_t disc, int_cl_t *ramified,
   int_cl_t *ord, cm_form_t *gen)
   /* Given a negative discriminant, compute a normal series of its ideal
      class group and return the length of the series. gen[i] generates a
      cyclic subgroup of prime order ord[i] of the classgroup divided by the
      group generated by gen[0],...,gen[i-1]. While the class group is not
      necessarily the direct product of the cyclic groups generated by
      gen[i], every class group element can be written uniquely as a
      product of gen[i]^e_i with 0<=e_i<ord[i]. So the result of this
      function can be used just like the Smith normal form of the class
      group to obtain all its elements, but it becomes useful when
      computing the class field as a tower of extensions.
      The array ord is sorted in decreasing order, which appears to be the
      best choice for the class field application.
      ord and gen must have been initialised with sufficient space to hold
      the result.
      If different from NULL, the 0-terminated array ramified [0], ...,
      ramified [n-1], ramified [n] = 0 must contain n>=2 ramified primes;
      then the function computes a normal series not for the class group,
      but for its quotient by the 2-elementary subgroup of order 2^(n-1)
      generated by the forms above p [0] * p [1], ..., p [0] * p [n-1]. */

{
   uint_cl_t factors [17];
   unsigned int exponents [17];
   int_cl_t o [64];
   cm_form_t g [64];
   uint_cl_t h, p;
   int l, p_l, length, i;

   if (ramified == NULL)
      l = cm_pari_classgroup (disc, o, g);
   else
      l = cm_pari_classgroup_2quotient (disc, ramified, o, g);

   if (l == 0) {
      length = 0;
   }
   else {
      h = 1;
      for (i = 0; i < l; i++)
         h *= o [i];

      cm_nt_factor (h, factors, exponents);
      for (p_l = 1; factors [p_l] != 0; p_l++);

      length = 0;
      /* Loop over all primes p dividing h. */
      for (p_l--; p_l >= 0; p_l--) {
         p = factors [p_l];
         /* Taking care of prime power factors requires a second loop. */
         while (h % p == 0) {
            /* From the first to the last cyclic factor of the classgroup,
               look for a cyclic subfactor of order p. */
            for (i = 0; i < l; i++)
               if (o [i] % p == 0) {
                  /* Record the cyclic subfactor. */
                  cm_classgroup_pow (&(gen [length]), g [i], o [i] / p, disc);
                  ord [length] = p;
                  length++;
                  /* Quotient by the cyclic subfactor. */
                  o [i] /= p;
                  h /= p;
               }
         }
      }
   }

   return length;
}

/*****************************************************************************/
/*****************************************************************************/
