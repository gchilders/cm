/*

nt.c - number theoretic helper functions

Copyright (C) 2009 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include "cm_common-impl.h"

static void nt_elliptic_curve_double (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t a, mpz_t p);
static void nt_elliptic_curve_add (mpz_t P_x,  mpz_t P_y, bool *P_infty,
   mpz_t Q_x, mpz_t Q_y, const bool Q_infty, mpz_t a, mpz_t p);

/*****************************************************************************/

long int cm_nt_mod (const long int a, const long int p)
   /* returns a representative of a % p in [0, p-1[ */

{
   long int res = a % p;
   return (res < 0 ? res + p : res);
}

/*****************************************************************************/

long int cm_nt_sqrt (const unsigned long int n)
   /* returns the positive root of n if it is a perfect square, -1 otherwise */

{
   double   root;
   unsigned long int s;

   root = sqrt ((double) n);
   s = (unsigned long int) (root + 0.5);
   if (s * s == n)
      return (long int) s;
   else
      return -1;
}

/*****************************************************************************/

long int cm_nt_gcd (long int a, long int b)
   /* returns the positive gcd of a and b */

{
   long int r, a_local, b_local;

   if (a == 0)
      return b;
   else if (b == 0)
      return a;
   else
   {
      if (a < 0)
         a_local = -a;
      else
         a_local = a;
      if (b < 0)
         b_local = -b;
      else
         b_local = b;
      r = a_local % b_local;
      while (r > 0)
      {
         a_local = b_local;
         b_local = r;
         r = a_local % b_local;
      }
      return b_local;
   }
}

/*****************************************************************************/

long int cm_nt_gcdext (long int *u, long int *v, long int a, long int b)
   /* returns the positive gcd d of a and b and u, v such that u a + v b = d */

{
   long int r0, r1, r2, u0, u1, u2, v0, v1, v2, q;
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

   *u = sgn_a * u0;
   *v = sgn_b * v0;

   return r0;
}

/*****************************************************************************/

long int cm_nt_invert (long int a, const long int p)
   /* returns a^(-1) mod p */

{
   long int r, a_local, b_local, u0, u1, u2;

   if (a == 0)
   {
      printf ("*** nt_invert called with 0\n");
      exit (1);
   }
   else if (a < 0)
      /* we assume that a is reduced mod p by the C function */
      a += p;

   a_local = a;
   b_local = p;
   u0 = 1;
   u1 = 0;
   /* loop invariant: u0 * a + something * p = a_local */
   /* loop invariant: u1 * a + something * p = b_local */
   while (b_local > 0)
   {
      u2 = u0 - (a_local / b_local) * u1;
      r = a_local % b_local;
      u0 = u1;
      u1 = u2;
      a_local = b_local;
      b_local = r;
   }
   return u0;
}

/*****************************************************************************/

int cm_nt_kronecker (long int a, long int b)
   /* computes the Kronecker symbol (a / b) following Algorithm 1.4.12       */
   /* in Cohen (by the binary algorithm)                                     */
{
   const int tab [8] = {0, 1, 0, -1, 0, -1, 0, 1};
   int k, r;

   /* step 1 */
   if (b == 0)
   {
      if (a == 1 || a == -1)
         return 1;
      else
         return 0;
   }

   /* step 2*/
   if (a % 2 == 0 && b % 2 == 0)
      return 0;

   while (b % 4 == 0)
      b >>= 2;
   if (b % 2 == 0)
   {
      b >>= 1;
      k = tab [a & 7];
   }
   else
      k = 1;

   if (b < 0)
   {
      b = -b;
      if (a < 0)
         k = -k;
   }

   /* step 3 and added test; here, b is already odd */
   if (a < 0)
   {
      a = -a;
      if (b & 2)
         k = -k;
   }
   a %= b;

   /* steps 4 to 6 */
   while (a != 0)
   {
      while (a % 4 == 0)
         a >>= 2;
      if (a % 2 == 0)
      {
         a >>= 1;
         k *= tab [b & 7];
      }
      if (b > a)
      {
         r = b - a;
         if (a & b & 2)
            k = -k;
         b = a;
         a = r;
      }
      else
         a -= b;
   }

   if (b > 1)
      return 0;
   else
      return k;
}

/*****************************************************************************/

int cm_nt_is_prime (mpz_t a)

{
   return (mpz_probab_prime_p (a, 10) > 0);
}

/*****************************************************************************/

int cm_nt_is_prime_l (const unsigned long int prime)

{
   int res;
   mpz_t a;

   mpz_init_set_ui (a, prime);
   res = cm_nt_is_prime (a);

   mpz_clear (a);

   return res;
}

/*****************************************************************************/

unsigned long int cm_nt_next_prime (const unsigned long int n)
   /* returns the prime following n */

{
   const unsigned long int P [] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
      31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
      107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
      181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257,
      263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
      349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
      433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509,
      521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
      613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
      701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
      809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883,
      887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991,
      997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063,
      1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
      1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237,
      1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319,
      1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433,
      1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
      1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597,
      1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669,
      1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777,
      1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873,
      1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979,
      1987, 1993, 1997, 1999 };
   const unsigned int size = 303; /* number of entries in P */

   if (n < P [size - 1]) {
      /* search in prime list; loop invariants:
         left <= right
         P [0], ..., P [left - 1] <= n
         P [right], ..., P [size - 1] > n      */
      unsigned int left = 0, right = size - 1, middle;
      while (left != right) {
         middle = (left + right) / 2; /* < right */
         if (P [middle] > n)
            right = middle; /* becomes smaller */
         else
            left = middle + 1; /* becomes bigger, remains <= right */
      }
      return P [left];
   }
   else {
      unsigned long int res;
      mpz_t a;
      if (n & 1)
         mpz_init_set_ui (a, n+2);
      else
         mpz_init_set_ui (a, n+1);
      while (!cm_nt_is_prime (a))
         mpz_add_ui (a, a, 2ul);
      res = mpz_get_ui (a);
      mpz_clear (a);
      return res;
   }
}

/*****************************************************************************/

void cm_nt_factor (long int d, unsigned long int *factors,
   unsigned int *exponents)
   /* factors the absolute value of d by trial division. The prime factors   */
   /* are stored in "factors", their multiplicities in "exponents", which    */
   /* must provide sufficiently much space. The list of prime factors is     */
   /* terminated by 0, so that 12 entries suffice for a number of 32 bits,   */
   /* and 17 entries for a number of 64 bits.                                */

{
   unsigned long int no, trial, trial2;
   int i, j;

   if (d < 0)
      no = -d;
   else
      no = d;

   i = 0;
   j = 0;
   trial = 0;
   trial2 = 0;
   while (trial2 <= no)
   {
      if (trial == 0)
      {
         trial = 2;
         trial2 = 4;
      }
      else if (trial == 2)
      {
         trial = 3;
         trial2 = 9;
      }
      else
      {
         trial += 2;
         trial2 += 4 * (trial - 1);
      }
      if (no % trial == 0)
      {
         factors [j] = trial;
         exponents [j] = 1;
         no /= trial;
         while (no % trial == 0)
         {
            no /= trial;
            exponents [j]++;
         }
         j++;
      }
      i++;
   }
   if (no != 1)
   {
     factors [j] = no;
     exponents [j] = 1;
     j++;
   }
   factors [j] = 0;
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_z (mpz_t root, mpz_t a, mpz_t p, mpz_t q)
   /* computes a square root of a modulo q=p^n, first modulo p by the        */
   /* Tonelli-Shanks algorithm, see Cohen, Algorithm 1.5., then doing a      */
   /* Hensel lift                                                            */

{
   mpz_t             pm1, z;
   mpz_t             a_local, y, x, b, tmp, tmp2, p_i;
   unsigned int      e;
   unsigned long int r, m;

   mpz_init (a_local);
   mpz_tdiv_r (a_local, a, q);
   if (mpz_cmp_ui (a_local, 0ul) == 0) {
      mpz_set_ui (root, 0ul);
      return;
   }

   mpz_init (pm1);
   mpz_init (z);
   mpz_init (y);
   mpz_init (x);
   mpz_init (b);
   mpz_init (tmp);
   mpz_init (tmp2);
   mpz_init (p_i);

   mpz_sub_ui (pm1, p, 1ul);
   e = 0;
   while (mpz_divisible_2exp_p (pm1, 1ul) != 0) {
      mpz_tdiv_q_2exp (pm1, pm1, 1ul);
      e++;
   }
   if (e > 1) {
      /* find generator of the 2-Sylow group */
      for (mpz_set_ui (z, 2ul); mpz_legendre (z, p) != -1;
               mpz_add_ui (z, z, 1ul));
      mpz_powm (z, z, pm1, p);
   }

   if (e == 1) /* p=3 (mod 8) */ {
      mpz_add_ui (tmp, p, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 2ul);
      mpz_powm (x, a_local, tmp, p);
   }
   else {
      /* initialisation */
      mpz_set (y, z);
      r = e;
      mpz_sub_ui (tmp, pm1, 1ul);
      mpz_tdiv_q_2exp (tmp, tmp, 1ul);
      mpz_powm (x, a, tmp, p);
      mpz_powm_ui (b, x, 2ul, p);
      mpz_mul (b, b, a_local);
      mpz_mod (b, b, p);
      mpz_mul (x, x, a_local);
      mpz_mod (x, x, p);
      while (mpz_cmp_ui (b, 1ul) != 0) {
         /* find exponent */
         m = 1;
         mpz_powm_ui (tmp, b, 2ul, p);
         while (mpz_cmp_ui (tmp, 1ul) != 0) {
            m++;
            mpz_powm_ui (tmp, tmp, 2ul, p);
         }
         if (m == r) {
            printf ("*** mpz_tonelli called with a = ");
            mpz_out_str (stdout, 10, a);
            printf (" and p = ");
            mpz_out_str (stdout, 10, p);
            printf (",\nbut a is not a square modulo p!\n");
            exit (1);
         }
         r = 1 << (r - m - 1);
         mpz_powm_ui (tmp, y, r, p);
         mpz_powm_ui (y, tmp, 2ul, p);
         r = m;
         mpz_mul (x, x, tmp);
         mpz_mod (x, x, p);
         mpz_mul (b, b, y);
         mpz_mod (b, b, p);
      }
   }

   mpz_set (p_i, p);

   while (mpz_cmp (p_i, q) < 0) {
      /* The root x is known modulo p_i, refine by Hensel for p_i^2. */
      mpz_mod (x, x, p_i);
      mpz_pow_ui (tmp, x, 2ul);
      mpz_sub (tmp, tmp, a_local);
      mpz_divexact (tmp, tmp, p_i);

      mpz_mul_2exp (tmp2, x, 1ul);
      mpz_invert (tmp2, tmp2, p_i);
      mpz_mul (tmp, tmp, tmp2);
      mpz_mod (tmp, tmp, p_i);

      mpz_mul (tmp, tmp, p_i);
      mpz_sub (x, x, tmp);

      mpz_pow_ui (p_i, p_i, 2ul);
   }

   mpz_mod (root, x, q);

   mpz_clear (a_local);
   mpz_clear (pm1);
   mpz_clear (z);
   mpz_clear (y);
   mpz_clear (x);
   mpz_clear (b);
   mpz_clear (tmp);
   mpz_clear (tmp2);
   mpz_clear (p_i);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli (mpz_t root, const long int a, mpz_t p, mpz_t q)
   /* computes a square root of a modulo q=p^n */

{
   mpz_t tmp_a;

   mpz_init_set_si (tmp_a, a);
   cm_nt_mpz_tonelli_z (root, tmp_a, p, q);
   mpz_clear (tmp_a);
}

/*****************************************************************************/
/*****************************************************************************/

static void nt_elliptic_curve_double (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t a, mpz_t p)
   /* replaces P by 2P on the elliptic curve given by a and the prime field  */
   /* of characteristic p; b is implicit since it is assumed that the point  */
   /* lies on the curve.                                                     */
   /* P_x and P_y are the X- and Y-coordinates of the point, respectively,   */
   /* P_infty is true if the point is the point at infinity. In this case,   */
   /* the values of P_x and P_y are not defined.                             */
   /* P_x, P_y and P_infty are modified.                                     */

{
   mpz_t factor, x_local, tmp, tmp2;

   mpz_init (factor);
   mpz_init (x_local);
   mpz_init (tmp);
   mpz_init (tmp2);

   if (!(*P_infty))
   /* there is nothing to do when P is already infinity */
   {
      if (mpz_cmp_ui (P_y, 0ul) == 0)
         *P_infty = true;
      else
      {
         /* factor = (3*P_x^2 + a) / 2*P_y */
         mpz_pow_ui (tmp, P_x, 2ul);
         mpz_mod (tmp, tmp, p);
         mpz_mul_ui (factor, tmp, 3ul);
         mpz_add (factor, factor, a);
         mpz_mul_2exp (tmp, P_y, 1ul);
         mpz_invert (tmp2, tmp, p);
         mpz_mul (tmp, factor, tmp2);
         mpz_mod (factor, tmp, p);
         /* P_x = factor^2 - 2 P_x; save P_x */
         mpz_set (x_local, P_x);
         mpz_pow_ui (P_x, factor, 2ul);
         mpz_submul_ui (P_x, x_local, 2ul);
         mpz_mod (P_x, P_x, p);
         /* P_y = factor * (old x - new x) - P_y */
         mpz_sub (tmp, x_local, P_x);
         mpz_mul (tmp2, factor, tmp);
         mpz_sub (P_y, tmp2, P_y);
         mpz_mod (P_y, P_y, p);
      }
   }

   mpz_clear (factor);
   mpz_clear (x_local);
   mpz_clear (tmp);
   mpz_clear (tmp2);
}

/*****************************************************************************/

static void nt_elliptic_curve_add (mpz_t P_x,  mpz_t P_y, bool *P_infty,
   mpz_t Q_x, mpz_t Q_y, const bool Q_infty, mpz_t a, mpz_t p)
   /* replaces P by P+Q on the elliptic curve given by a and the implicit b  */
   /* over the prime field of characteristic p                               */
   /* P_x, P_y and P_infty are modified.                                     */

{
   mpz_t factor, x_local, tmp, tmp2;

   mpz_init (factor);
   mpz_init (x_local);
   mpz_init (tmp);
   mpz_init (tmp2);

   if (!Q_infty)
   {
      if (*P_infty)
      {
         mpz_set (P_x, Q_x);
         mpz_set (P_y, Q_y);
         *P_infty = false;
      }
      else if (mpz_cmp (P_x, Q_x) == 0)
      {
         if (mpz_cmp (P_y, Q_y) == 0)
            nt_elliptic_curve_double (P_x, P_y, P_infty, a, p);
         else
            *P_infty = true;
      }
      else
      {
         /* factor = Delta y / Delta x */
         mpz_sub (tmp, Q_x, P_x);
         mpz_invert (tmp2, tmp, p);
         mpz_sub (tmp, Q_y, P_y);
         mpz_mul (factor, tmp, tmp2);
         mpz_mod (factor, factor, p);
         /* P_x = factor^2 - (P_x + Q_x), keep a copy of P_x */
         mpz_set (x_local, P_x);
         mpz_pow_ui (tmp, factor, 2ul);
         mpz_add (tmp2, P_x, Q_x);
         mpz_sub (P_x, tmp, tmp2);
         mpz_mod (P_x, P_x, p);
         /* P_y = factor * (old P_x - new x) - P_y */
         mpz_sub (tmp, x_local, P_x);
         mpz_mul (tmp2, factor, tmp);
         mpz_sub (P_y, tmp2, P_y);
         mpz_mod (P_y, P_y, p);
      }
   }

   mpz_clear (factor);
   mpz_clear (x_local);
   mpz_clear (tmp);
   mpz_clear (tmp2);
}

/*****************************************************************************/

void cm_nt_elliptic_curve_multiply (mpz_t P_x, mpz_t P_y, bool *P_infty,
   mpz_t m, mpz_t a, mpz_t p)
   /* replaces P by mP on the elliptic curve given by a and the implicit b   */
   /* over the prime field of characteristic p                               */
   /* P_x, P_y and P_infty are modified.                                     */

{
   mpz_t Q_x, Q_y, m_local, m_new;
   bool  Q_infty = true;
      /* m P is stored in Q; we use a right to left binary exponentiation    */
      /* scheme with the invariant Q + m P = const.                          */

   mpz_init (Q_x);
   mpz_init (Q_y);
   mpz_init (m_local);
   mpz_init (m_new);
   *P_infty = false;

   mpz_set (m_local, m);
   if (mpz_cmp_ui (m_local, 0ul) == 0)
      *P_infty = true;
   else
   {
      if (mpz_cmp_ui (m_local, 0ul) < 0)
      {
         mpz_neg (m_local, m_local);
         mpz_neg (P_y, P_y);
      }
      while (mpz_cmp_ui (m_local, 0ul) > 0)
      {
         while (mpz_tdiv_q_ui (m_new, m_local, 2ul) == 0)
         {
            mpz_set (m_local, m_new);
            nt_elliptic_curve_double (P_x, P_y, P_infty, a, p);
         }
         mpz_sub_ui (m_local, m_local, 1ul);
         nt_elliptic_curve_add (Q_x, Q_y, &Q_infty, P_x, P_y, *P_infty, a, p);
      }
      /* copy the result */
      mpz_set (P_x, Q_x);
      mpz_set (P_y, Q_y);
      *P_infty = Q_infty;
   }

   mpz_clear (Q_x);
   mpz_clear (Q_y);
   mpz_clear (m_local);
   mpz_clear (m_new);
}

/*****************************************************************************/

void cm_nt_elliptic_curve_random (mpz_t P_x, mpz_t P_y,
   mpz_t cofactor, mpz_t a, mpz_t b, mpz_t p)
   /* creates a point on the elliptic curve given by a and b over the prime  */
   /* field of characteristic p and multiplies it by the cofactor until      */
   /* the result is different from infinity. If the curve order is cofactor  */
   /* times a prime, this results in a point of order this prime.            */
   /* The point is not really random, since successive X-coordinates from    */
   /* 1 on are tested.                                                       */
   /* P_x and P_y are modified.                                              */

{
   mpz_t  tmp;
   long unsigned int P_x_long = 0;
   bool P_infty = true;

   mpz_init (tmp);
   while (P_infty)
   {
      P_x_long++;
      /* P_y = P_x^3 + a P_x + b */
      mpz_mul_ui (P_y, a, P_x_long);
      mpz_add (P_y, P_y, b);
      mpz_add_ui (P_y, P_y, P_x_long * P_x_long * P_x_long);
      mpz_mod (P_y, P_y, p);
      /* try to compute the square root of P_y */
      if (mpz_jacobi (P_y, p) != -1)
      {
         mpz_set_ui (P_x, P_x_long);
         cm_nt_mpz_tonelli_z (P_y, P_y, p, p);
         /* get rid of the cofactor */
         P_infty = false;
         cm_nt_elliptic_curve_multiply (P_x, P_y, &P_infty, cofactor, a, p);
      }
   }
   mpz_clear (tmp);
}

/*****************************************************************************/
/*****************************************************************************/

bool cm_nt_mpfr_get_z (mpz_t out, mpfr_t in)
   /* tries to round the real to an integer value                            */
   /* The return value reflects the success of the operation.                */

{
   mpfr_t rounded, diff;
   bool ok;
   mp_exp_t exp;

   mpfr_init2 (rounded, mpfr_get_prec (in));
   mpfr_init2 (diff, mpfr_get_prec (in));

   mpfr_round (rounded, in);
   mpfr_sub (diff, in, rounded, GMP_RNDN);
   if (mpfr_sgn (diff) == 0 || (- diff->_mpfr_exp) >= 10) {
      exp = mpfr_get_z_exp (out, rounded);
      if (exp > 0)
         ok = false;
      else if (exp < 0) {
         mpz_tdiv_q_2exp (out, out, (unsigned long int) (-exp));
         ok = true;
      }
      else
         ok = true;
   }
   else
      ok = false;

   mpfr_clear (rounded);
   mpfr_clear (diff);

   return ok;
}

/*****************************************************************************/
/*****************************************************************************/
