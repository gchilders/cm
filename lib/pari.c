/*

pari.c - functions using pari; for factoring polynomials and for computing
generators of class groups

Copyright (C) 2010, 2015, 2018, 2021 Andreas Enge

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

static GEN mpz_get_Z (mpz_srcptr z);
static void Z_get_mpz (mpz_ptr z, GEN x);
static GEN mpzx_get_FpX (mpzx_srcptr f, mpz_srcptr p);

/*****************************************************************************/
/*                                                                           */
/* Conversion functions between PARI, GMP and internal types.                */
/*                                                                           */
/*****************************************************************************/

static GEN mpz_get_Z (mpz_srcptr z)
   /* returns the GEN of type t_INT corresponding to z */

{
   const long l = z->_mp_size;
   const long lz = labs (l) + 2;
   int i;
   GEN x = cgeti (lz);

   x [1] = evalsigne ((l > 0 ? 1 : (l < 0 ? -1 : 0))) | evallgefint (lz);
   for (i = 0; i < labs (l); i++)
      *int_W (x, i) = (z->_mp_d) [i];

   return x;
}

/*****************************************************************************/

static void Z_get_mpz (mpz_ptr z, GEN x)
   /* returns via z the gmp integer corresponding to the t_INT x */

{
   const long l = lgefint (x) - 2;
   int i;

   _mpz_realloc (z, l);
   z->_mp_size = (signe (x) > 0 ? l : -l);
   for (i = 0; i < l; i++)
      (z->_mp_d) [i] = *int_W (x, i);
}

/*****************************************************************************/

static GEN icl_get_Z (int_cl_t z)
   /* returns the GEN of type t_INT corresponding to z; for the time being,
      we use a quick and dirty implementation int_cl_t -> mpz_t -> GEN,
      while it should be easier to move the 64 bits around... */

{
   mpz_t zz;
   GEN x;

   mpz_init (zz);
   cm_classgroup_mpz_set_icl (zz, z);
   x = mpz_get_Z (zz);
   mpz_clear (zz);

   return x;
}

/*****************************************************************************/

static int_cl_t Z_get_icl (GEN x)
   /* returns the int_cl_t correspondong to x, again using a quick and
      dirty implementation */

{
   mpz_t zz;
   int_cl_t i;

   mpz_init (zz);
   Z_get_mpz (zz, x);
   i = cm_classgroup_mpz_get_icl (zz);
   mpz_clear (zz);

   return i;
}

/*****************************************************************************/

static GEN mpzx_get_FpX (mpzx_srcptr f, mpz_srcptr p)
   /* Return a GEN of type t_POL over t_INT corresponding to the polynomial
      f, reduced modulo p. */

{
   int i;
   GEN res;
   mpz_t tmp;

   mpz_init (tmp);

   res = cgetg (f->deg + 3, t_POL);
   setvarn (res, 0);
   if (f->deg == -1)
      setsigne (res, 0);
   else {
      setsigne (res, 1);
      for (i = 0; i <= f->deg; i++) {
         mpz_mod (tmp, f->coeff [i], p);
         gel (res, i+2) = mpz_get_Z (tmp);
      }
   }

   mpz_clear (tmp);

   return res;
}


/*****************************************************************************/
/*                                                                           */
/* Exported functions.                                                       */
/*                                                                           */
/*****************************************************************************/

void cm_pari_init ()
{
   pari_init_opts (2000000, 0, INIT_JMPm | INIT_DFTm);
      /* Do not capture SIGSEGV. */
   paristack_setsize (2000000, 1000000000);
}

/*****************************************************************************/

void cm_pari_clear ()

{
   pari_close ();
}

/*****************************************************************************/

bool cm_pari_eval_int (mpz_ptr n, char *e)
   /* If the PARI expression e evaluates to a t_INT, set n to this integer
      and return true; otherwise keep n unchanged and return false. */
{
   pari_sp av;
   GEN z;
   bool ok = true;

   av = avma;

   z = gp_read_str (e);
   if (typ (z) == t_INT)
      Z_get_mpz (n, z);
   else
      ok = false;

   avma = av;

   return ok;
}


/*****************************************************************************/
/*                                                                           */
/* Functions for finding roots of polynomials.                               */
/*                                                                           */
/*****************************************************************************/

void cm_pari_oneroot (mpz_ptr root, mpzx_srcptr f, mpz_srcptr p, bool verbose)
   /* Find a root of the polynomial f over the prime field of
      characteristic p, assuming that f splits completely, and return it
      in the variable of the same name. */

{
   pari_sp av;
   GEN fp, pp, rootp;
   cm_timer_t clock;

   av = avma;

   cm_timer_start (clock);
   if (verbose)
      printf ("--- Root finding in degree %i\n", f->deg);

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);

   rootp = FpX_oneroot_split (fp, pp);
   Z_get_mpz (root, rootp);

   cm_timer_stop (clock);
   if (verbose)
      printf ("-- Time for root: %.1f\n", cm_timer_get (clock));

   avma = av;
}

/*****************************************************************************/

mpz_t* cm_pari_find_roots (int *no, mpzx_srcptr f, mpz_srcptr p)
   /* Computes all the roots (without multiplicities) of the polynomial f
      modulo p. The number of found roots is returned in no. */

{
   pari_sp av;
   mpz_t *res;
   GEN fp, pp, rootsp;
   int i;

   av = avma;

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);
   rootsp = FpX_roots (fp, pp);
   *no = lg (rootsp) - 1;
   res = (mpz_t*) malloc ((*no) * sizeof (mpz_t));
   for (i = 0; i < *no; i++) {
      mpz_init (res [i]);
      Z_get_mpz (res [i], gel (rootsp, i+1));
   }

   avma = av;

   return res;
}

/*****************************************************************************/
/*                                                                           */
/* Functions for computing class groups.                                     */
/*                                                                           */
/*****************************************************************************/

int cm_pari_classgroup (int_cl_t disc, int_cl_t *ord, cm_form_t *gen)
   /* Given a negative discriminant, compute its ideal class group as a
      product of cyclic groups with their orders and generators. The orders
      are returned in ord and the generators in gen; their number is the
      return value. ord and gen must have been initialised with sufficient
      space to hold the result. */

{
   pari_sp av;
   int length;
   GEN d, cl, orders, gens, qfb;
   int i;

   av = avma;

   d = icl_get_Z (disc);
   cl = quadclassunit0 (d, 0, NULL, 0);
   orders = gel (cl, 2);
   gens = gel (cl, 3);
   length = glength (orders);
   for (i = 0; i < length; i++) {
      ord [i] = Z_get_icl (gel (orders, i+1));
      qfb = gel (gens, i+1);
      gen [i].a = Z_get_icl (gel (qfb, 1));
      gen [i].b = Z_get_icl (gel (qfb, 2));
   }

   avma = av;

   return length;
}

/*****************************************************************************/

int cm_pari_classgroup_2quotient (int_cl_t disc, const int *p,
   int_cl_t *ord, cm_form_t *gen)
   /* Given a negative discriminant disc and a 0-terminated list of primes
      p [0], p [1], ..., p [n-1] dividing the fundamental part of disc
      with n>=2, compute the quotient of the ideal class group by the
      subgroup of order 2^(n-1) generated by the forms of norm
      p [0] * p [1], ..., p [0] * p [n-1], as a product of cyclic groups
      with their orders and generators. The orders are returned in ord and
      the generators in gen; their number is the return value. ord and gen
      must have been initialised with sufficient space to hold the result.
      The function could be combined with the previous one, but is kept
      apart for the sake of the clarity of the previous function. */

{
   pari_sp av;
   int length;
   GEN d, cl, orders, gens, qfb;
   int length2, size2;
   GEN gens2, group2, dlog2, halforder;
   GEN M, Uinv, ordersq, gensq;
   int lengthq;
   cm_form_t p0pi;
   int i, j;

   av = avma;

   /* Compute the form class group. */
   d = icl_get_Z (disc);
   cl = quadclassunit0 (d, 0, NULL, 0);
   orders = gel (cl, 2);
   gens = gel (cl, 3);
   length = glength (orders);

   /* Compute the size of the 2-elementary subgroup and its generators
      as quadratic forms. */
   for (length2 = 0;
        length2 < length && !mpodd (gel (orders, length2 + 1));
        length2++);
   gens2 = cgetg (length2 + 1, t_VEC);
   for (i = 1; i <= length2; i++)
      gel (gens2, i) = gpowgs (gel (gens, i), itou (gel (orders, i)) / 2);

   /* Enumerate the subgroup elements and their discrete logarithms. */
   size2 = 1 << length2;
   group2 = cgetg (size2 + 1, t_VEC);
   dlog2 = cgetg (size2 + 1, t_VEC);
   size2 = 1;
#if PARI_VERSION_CODE < PARI_VERSION (2, 14, 0)
   gel (group2, 1) = qfi_1 (gel (gens, 1));
#else
   gel (group2, 1) = qfb_1 (gel (gens, 1));
#endif
   gel (dlog2, 1) = zerocol (length);
   for (i = 1; i <= length2; i++) {
      halforder = shifti (gel (orders, i), -1);
      for (j = 1; j <= size2; j++) {
#if PARI_VERSION_CODE < PARI_VERSION (2, 14, 0)
         gel (group2, size2 + j)
            = qficomp (gel (group2, j), gel (gens2, i));
#else
         gel (group2, size2 + j)
            = qfbcomp (gel (group2, j), gel (gens2, i));
#endif
         gel (dlog2, size2 + j) = shallowcopy (gel (dlog2, j));
         gmael (dlog2, size2 + j, i) = halforder;
      }
      size2 <<= 1;
   }

   /* Create a matrix of relations for the quotient group. */
   M = zeromat (length, 0);
   if (disc % p [0] != 0) {
      printf ("***** Error: Calling cm_pari_classgroup_2quotient with "
              "non-ramified prime %i.\n", p [0]);
      exit (1);
   }
   for (i = 1; p [i] != 0; i++) {
      if (disc % p [i] != 0) {
         printf ("***** Error: Calling cm_pari_classgroup_2quotient with "
                 "non-ramified prime %i.\n", p [i]);
         exit (1);
      }
      /* Compute the form of norm p [0] * p [i]. */
      p0pi.a = p [0] * p [i];
      if (disc % 2 != 0)
         p0pi.b = p0pi.a;
      else if (disc % 8 == 0 || p0pi.a % 2 != 0)
         p0pi.b = 0;
      else
         p0pi.b = p0pi.a;
      cm_classgroup_reduce (&p0pi, disc);
#if PARI_VERSION_CODE < PARI_VERSION (2, 14, 0)
      qfb = qfi (icl_get_Z (p0pi.a), icl_get_Z (p0pi.b),
         icl_get_Z (cm_classgroup_compute_c (p0pi.a, p0pi.b, disc)));
#else
      qfb = mkqfb (icl_get_Z (p0pi.a), icl_get_Z (p0pi.b),
         icl_get_Z (cm_classgroup_compute_c (p0pi.a, p0pi.b, disc)),
         icl_get_Z (disc));
#endif
      /* Look up its discrete logarithm. */
      for (j = 1; !gequal (qfb, gel (group2, j)); j++);
      M = shallowconcat (M, gel (dlog2, j));
   }

   /* Compute the quotient group. */
   M = hnfmodid (M, orders);
   ordersq = ZM_snf_group (M, NULL, &Uinv);
   lengthq = glength (ordersq);
   gensq = cgetg (lengthq + 1, t_VEC);
   for (i = 1; i <= lengthq; i++)
      gel (gensq, i) = factorback2 (gens, gel (Uinv, i));

   for (i = 0; i < lengthq; i++) {
      ord [i] = Z_get_icl (gel (ordersq, i+1));
      qfb = gel (gensq, i+1);
      gen [i].a = Z_get_icl (gel (qfb, 1));
      gen [i].b = Z_get_icl (gel (qfb, 2));
   }

   avma = av;

   return lengthq;
}

/*****************************************************************************/
/*                                                                           */
/* Functions used for ECPP                                                   */
/*                                                                           */
/*****************************************************************************/

void cm_pari_prime_product (mpz_ptr prim, unsigned long int a,
   unsigned long int b)
   /* Return in prim the product of all primes p with a < p <= b. */
{
   pari_sp av;
   GEN res;

   av = avma;
   res = zv_prod_Z (primes_interval_zv (a+1, b));
   Z_get_mpz (prim, res);
   avma = av;
}

/*****************************************************************************/

bool cm_pari_cornacchia (mpz_ptr t, mpz_ptr v, mpz_srcptr p,
   mpz_srcptr root, const int_cl_t d)
   /* The function has the same interface as cm_nt_mpz_cornacchia. As PARI
      uses a half gcd, it should be asymptotically faster. */
{
   pari_sp av;
   GEN px, py;
   long res;

   av = avma;
   res = cornacchia2_sqrt (icl_get_Z (-d), mpz_get_Z (p), mpz_get_Z (root),
      &px, &py);
   if (res == 1) {
      Z_get_mpz (t, px);
      if (v != NULL)
         Z_get_mpz (v, py);
   }
   avma = av;

   return (res == 1);
}

/*****************************************************************************/

bool cm_pari_ecpp_check (mpz_t **cert, int depth)
   /* Given a complete ECPP certificate (after step 2) of length depth in
      cert, use PARI to check it and return its validity. */
{
   pari_sp av;
   bool res;
   GEN c, ci;
   int i, j;

   av = avma;

   c = cgetg (depth + 1, t_VEC);
   for (i = 0; i < depth; i++) {
      ci = cgetg (6, t_VEC);
      gel (c, i+1) = ci;
      for (j = 0; j < 4; j++)
         gel (ci, j+1) = mpz_get_Z (cert [i][j]);
      gel (ci, 5) = mkvec2 (mpz_get_Z (cert [i][4]),
                            mpz_get_Z (cert [i][5]));
   }

   res = (primecertisvalid (c) == 1);

   avma = av;

   return res;
}

/*****************************************************************************/
