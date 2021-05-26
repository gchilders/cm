/*

nt.c - number theoretic helper functions

Copyright (C) 2009, 2010, 2015, 2021 Andreas Enge

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

int cm_nt_kronecker (int_cl_t a, int_cl_t b)
   /* Compute the Kronecker symbol (a / b) following Algorithm 1.4.12
      in Cohen93 (by the binary algorithm). */
{
   const int tab [8] = {0, 1, 0, -1, 0, -1, 0, 1};
   int k;
   int_cl_t r;

   /* step 1 */
   if (b == 0) {
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
   if (b % 2 == 0) {
      b >>= 1;
      k = tab [cm_classgroup_mod (a, (uint_cl_t) 8)];
   }
   else
      k = 1;

   if (b < 0) {
      b = -b;
      if (a < 0)
         k = -k;
   }

   /* step 3 */
   a = cm_classgroup_mod (a, b);

   /* steps 4 to 6 */
   while (a != 0) {
      while (a % 4 == 0)
         a >>= 2;
      if (a % 2 == 0) {
         a >>= 1;
         k *= tab [cm_classgroup_mod (b, (uint_cl_t) 8)];
      }
      if (b > a) {
         r = b - a;
         if (   cm_classgroup_mod (a, (uint_cl_t) 4) == 3
             && cm_classgroup_mod (b, (uint_cl_t) 4) == 3)
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

void cm_nt_factor (uint_cl_t d, uint_cl_t *factors, unsigned int *exponents)
   /* Factor d by trial division. The prime factors are stored in factors,
      their multiplicities in exponents, which must provide sufficiently
      much space. The list of prime factors is terminated by 0, so that
      17 entries suffice for the 64 bit type. */

{
   uint_cl_t p, p2;
   int i;

   i = 0;
   p = 2;
   p2 = 4;
   do {
      if (d % p == 0) {
         factors [i] = p;
         exponents [i] = 1;
         d /= p;
         while (d % p == 0) {
            exponents [i]++;
            d /= p;
         }
         i++;
      }
      /* We may wish to use cm_nt_next_prime, but its implementation
         is slow beyond the precomputed table. */
      if (p == 2) {
         p = 3;
         p2 = 9;
      }
      else {
         p2 += 4*(p+1);
         p += 2;
      }
   } while (p2 <= d);

   if (d != 1) {
     /* There is a prime factor left. */
     factors [i] = d;
     exponents [i] = 1;
     i++;
   }

   factors [i] = 0;
}

/*****************************************************************************/

void cm_nt_mpz_tonelli_z (mpz_ptr root, mpz_srcptr a, mpz_srcptr p)
   /* computes a square root of a modulo p by the Tonelli-Shanks algorithm,  */
   /* see Cohen93, Algorithm 1.5.                                            */

{
   mpz_t             pm1, z;
   mpz_t             a_local, y, x, b, tmp;
   unsigned int      e;
   unsigned long int r, m;

   mpz_init (a_local);
   mpz_tdiv_r (a_local, a, p);
   if (mpz_cmp_ui (a_local, 0ul) == 0) {
      mpz_set_ui (root, 0ul);
      mpz_clear (a_local);
      return;
   }

   mpz_init (pm1);
   mpz_init (z);
   mpz_init (y);
   mpz_init (x);
   mpz_init (b);
   mpz_init (tmp);

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
      mpz_powm (x, a_local, tmp, p);
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

   mpz_set (root, x);

   mpz_clear (a_local);
   mpz_clear (pm1);
   mpz_clear (z);
   mpz_clear (y);
   mpz_clear (x);
   mpz_clear (b);
   mpz_clear (tmp);
}

/*****************************************************************************/

void cm_nt_mpz_tonelli (mpz_ptr root, const long int a, mpz_srcptr p)
   /* computes a square root of a modulo p */

{
   mpz_t tmp_a;

   mpz_init_set_si (tmp_a, a);
   cm_nt_mpz_tonelli_z (root, tmp_a, p);
   mpz_clear (tmp_a);
}

/*****************************************************************************/

bool cm_nt_mpz_cornacchia (mpz_ptr t, mpz_ptr v, mpz_srcptr p,
   const int_cl_t d)
   /* Compute t such that 4*p = t^2-v^2*d for some v, where p is an odd
      prime and d is an imaginary-quadratic discriminant such that d is a
      square modulo p and |d|<4*p, using Algorithm 1.5.3 of Cohen93.
      The return value indicates whether such a t exists; if not, the
      value of t is not changed during the algorithm. If yes and v is not
      NULL, it is changed. */
{
   mpz_t r, a, b, l;
   bool ok;

   mpz_init (r);
   mpz_init (a);
   mpz_init (b);
   mpz_init (l);

   /* Step 3: Initialisation. */
   cm_nt_mpz_tonelli (b, d, p);
   if (mpz_divisible_2exp_p (b, 1)) {
      if (d % 4 != 0)
         mpz_sub (b, p, b);
   }
   else {
      if (d % 4 == 0)
         mpz_sub (b, p, b);
   }
   mpz_mul_2exp (l, p, 2);
   mpz_sqrt (l, l);
   mpz_mul_2exp (a, p, 1);

   /* Step 4: Euclidian algorithm. */
   while (mpz_cmp (b, l) > 0) {
      mpz_mod (r, a, b);
      mpz_set (a, b);
      mpz_set (b, r);
   }

   /* Step 5: Check correctness. */
   mpz_mul_2exp (a, p, 2);
   mpz_pow_ui (r, b, 2);
   mpz_sub (a, a, r); /* a = 4*p-b^2 */
   if (!(mpz_divisible_ui_p (a, -d)))
      ok = false;
   else {
      mpz_divexact_ui (a, a, -d);
      if (v == NULL) {
         if (!mpz_perfect_square_p (a))
            ok = false;
         else {
            ok = true;
            mpz_set (t, b);
         }
      }
      else {
         if (!mpz_root (l, a, 2))
            ok = false;
         else {
            ok = true;
            mpz_set (t, b);
            mpz_set (v, l);
         }
      }
   }

   mpz_clear (r);
   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (l);

   return ok;
}

/*****************************************************************************/
/*****************************************************************************/

bool cm_nt_fget_z (mpz_t out, ftype in)
   /* Tries to round the real value in to an integer value out. The return
      value reflects the success of the operation. */

{
   ftype rounded, diff;
   bool ok;
   mp_exp_t expo;

   finit (rounded, fget_prec (in));
   finit (diff, fget_prec (in));

   fround (rounded, in);
   fsub (diff, in, rounded);
   if (fsgn (diff) == 0 || (- fget_exp (diff)) >= 10) {
      expo = fget_z_exp (out, rounded);
      if (expo > 0)
         ok = false;
      else if (expo < 0) {
         mpz_tdiv_q_2exp (out, out, (unsigned long int) (-expo));
         ok = true;
      }
      else
         ok = true;
   }
   else
      ok = false;

   if (!ok) {
      printf ("***** Error in rounding:\n");
      fout_str (stdout, 10, 0ul, in);
      printf ("\n");
      fout_str (stdout, 10, 0ul, rounded);
      printf ("\n");
   }

   fclear (rounded);
   fclear (diff);

   return ok;
}

/*****************************************************************************/

bool cm_nt_cget_zz (mpz_ptr out1, mpz_ptr out2, ctype in, ctype omega)
   /* Tries to round the complex value to an imaginary-quadratic integer
      out1+out2*omega, where omega is the second element of an integral
      basis (or more generally, a non-real complex number). The return
      value reflects the success of the operation. */
{
   ftype tmp;
   bool ok;

   finit (tmp, cget_prec (in));

   fdiv (tmp, cimagref (in), cimagref (omega));
   ok = cm_nt_fget_z (out2, tmp);

   if (ok) {
      fmul_z (tmp, crealref (omega), out2);
      fsub (tmp, crealref (in), tmp);
      ok = cm_nt_fget_z (out1, tmp);
   }

   fclear (tmp);

   return ok;
}

/*****************************************************************************/
/*****************************************************************************/
