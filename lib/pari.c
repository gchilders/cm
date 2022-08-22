/*

pari.c - functions using pari; for factoring polynomials and for computing
generators of class groups

Copyright (C) 2010, 2015, 2018, 2021, 2022 Andreas Enge

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
static GEN icl_get_Z (int_cl_t z);
static int_cl_t Z_get_icl (GEN x);
static GEN mpzx_get_FpX (mpzx_srcptr f, mpz_srcptr p);
static void FpX_get_mpzx (mpzx_ptr f, GEN x);
static void mpzx_pow_modmod (mpzx_ptr g, mpzx_srcptr f, mpz_srcptr e,
   mpzx_srcptr m, mpz_srcptr p);
static void mpzx_gcd_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p);
static void mpzx_divexact_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p);
static int good_root_of_unity (mpz_ptr zeta, mpz_srcptr p, const int deg);
static void mpzx_oneroot_split_mod_rec (mpz_ptr root, mpzx_srcptr f,
   mpz_srcptr p, bool verbose, bool debug);

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
   /* Return a GEN of type t_POL over t_INT corresponding to the
      polynomial f, reduced modulo p. */

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

static void FpX_get_mpzx (mpzx_ptr f, GEN x)
   /* Return the lift of x from FpX to a polynomial over the integers
      in f. */
{
   int deg, i;

   deg = (signe (x) == 0 ? -1 : (int) poldegree (x, 0));
   mpzx_set_deg (f, deg);
   for (i = 0; i <= deg; i++)
      Z_get_mpz (f->coeff [i], gel (x, i+2));
}

/*****************************************************************************/
/*                                                                           */
/* Functions for mpzx modulo p relying on PARI.                              */
/*                                                                           */
/*****************************************************************************/

static void mpzx_pow_modmod (mpzx_ptr g, mpzx_srcptr f, mpz_srcptr e,
   mpzx_srcptr m, mpz_srcptr p)
   /* Compute g = f^e modulo m and p, assuming e >= 0. */
{
#ifdef HAVE_FLINT
   fmpz_t pp, ep;
   fmpz_mod_ctx_t ctx;
   fmpz_mod_poly_t fp, mp, gp;

   fmpz_init (pp);
   fmpz_set_mpz (pp, p);
   fmpz_init (ep);
   fmpz_mod_ctx_init (ctx, pp);
   fmpz_mod_poly_init (fp, ctx);
   fmpz_mod_poly_init (mp, ctx);
   fmpz_mod_poly_init (gp, ctx);

   fmpz_set_mpz (ep, e);
   fmpz_mod_poly_set_mpzx (fp, f, ctx);
   fmpz_mod_poly_set_mpzx (mp, m, ctx);

   fmpz_mod_poly_powmod_fmpz_binexp (gp, fp, ep, mp, ctx);

   mpzx_set_fmpz_mod_poly (g, gp, ctx);

   fmpz_clear (pp);
   fmpz_clear (ep);
   fmpz_mod_poly_clear (fp, ctx);
   fmpz_mod_poly_clear (mp, ctx);
   fmpz_mod_poly_clear (gp, ctx);
   fmpz_mod_ctx_clear (ctx);
#else
   GEN pp, ep, fp, mp, gp;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   ep = mpz_get_Z (e);
   fp = mpzx_get_FpX (f, p);
   mp = mpzx_get_FpX (m, p);

   gp = FpXQ_pow (fp, ep, mp, pp);

   FpX_get_mpzx (g, gp);

   avma = av;
#endif
}

/*****************************************************************************/

static void mpzx_gcd_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p)
   /* Compute h = gcd (f, g) modulo p. */
{
#ifdef HAVE_FLINT
   fmpz_t pp;
   fmpz_mod_ctx_t ctx;
   fmpz_mod_poly_t fp, gp, hp;

   fmpz_init (pp);
   fmpz_set_mpz (pp, p);
   fmpz_mod_ctx_init (ctx, pp);
   fmpz_mod_poly_init (fp, ctx);
   fmpz_mod_poly_init (gp, ctx);
   fmpz_mod_poly_init (hp, ctx);

   fmpz_mod_poly_set_mpzx (fp, f, ctx);
   fmpz_mod_poly_set_mpzx (gp, g, ctx);

   fmpz_mod_poly_gcd (hp, fp, gp, ctx);

   mpzx_set_fmpz_mod_poly (h, hp, ctx);

   fmpz_clear (pp);
   fmpz_mod_poly_clear (fp, ctx);
   fmpz_mod_poly_clear (gp, ctx);
   fmpz_mod_poly_clear (hp, ctx);
   fmpz_mod_ctx_clear (ctx);
#else
   GEN pp, fp, gp, hp;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);
   gp = mpzx_get_FpX (g, p);

   hp = FpX_gcd (fp, gp, pp);

   FpX_get_mpzx (h, hp);

   avma = av;
#endif
}

/*****************************************************************************/

static void mpzx_divexact_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p)
   /* Assuming that g divides f, compute the quotient in h. */
{
   GEN fp, gp, pp, hp;

   pari_sp av = avma;

   fp = mpzx_get_FpX (f, p);
   gp = mpzx_get_FpX (g, p);
   pp = mpz_get_Z (p);

   hp = FpX_div (fp, gp, pp);

   FpX_get_mpzx (h, hp);

   avma = av;
}

/*****************************************************************************/
/*                                                                           */
/* Various simple functions.                                                 */
/*                                                                           */
/*****************************************************************************/

void cm_pari_init ()
{
   pari_init_opts (1ul<<23, 0, INIT_JMPm | INIT_DFTm);
      /* Do not capture SIGSEGV. */
   paristack_setsize (1ul<<23, 1ul<<31);
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

char* cm_pari_sprintf_hfactor (int_cl_t d)
   /* Return a newly allocated string containing the factorisation of the
      class number of d. */
{
   char *h, *res;
   pari_sp av = avma;

   default0 ("output", "0");
   h = pari_sprintf ("%Ps", factorint (qfbclassno0 (icl_get_Z (d), 0), 0));
   res = (char *) malloc ((strlen (h) + 1) * sizeof (char));
   strcpy (res, h);

   avma = av;

   return res;
}

/*****************************************************************************/
/*                                                                           */
/* Functions for finding roots of polynomials.                               */
/*                                                                           */
/*****************************************************************************/

void cm_pari_oneroot (mpz_ptr root, mpzx_srcptr f, mpz_srcptr p)
   /* Find a root of the polynomial f over the prime field of
      characteristic p, assuming that f splits completely, and return it
      in the variable of the same name. */
{
   GEN fp, pp, rootp;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, p);

   rootp = FpX_oneroot_split (fp, pp);
   Z_get_mpz (root, rootp);

   avma = av;
}

/*****************************************************************************/

static int good_root_of_unity (mpz_ptr zeta, mpz_srcptr p, const int deg)
   /* Compute in zeta a root of unity in F_p with p odd that is suitable
      for finding a root of a totally split polynomial of degree deg > 1;
      its order n is returned. A good choice seems to be n close to deg;
      we decrement it until it divides p-1. */

{
   GEN pp, pm1, factn, e, base, zetap;
   int n;

   pari_sp av = avma;

   pp = mpz_get_Z (p);
   pm1 = subis (pp, 1ul);
   for (n = deg; !dvdiu (pm1, n); n--);

   factn = Z_factor (stoi (n));
   e = diviuexact (pm1, n);
   base = gen_1;
   do {
      base = addis (base, 1l);
      zetap = Fp_pow (base, e, pp);
   }
   while (!equaliu (Fp_order (zetap, factn, pp), n));

   Z_get_mpz (zeta, zetap);

   avma = av;

   return n;
}

/*****************************************************************************/

static void mpzx_oneroot_split_mod_rec (mpz_ptr root, mpzx_srcptr f,
   mpz_srcptr p, bool verbose, bool debug)
   /* Compute in root a root of the polynomial f over the prime field
      of characteristic p, assuming that f splits completely and that
      its coefficients are reduced modulo p. */
{
   int n, target, min, i;
   mpz_t zeta, e, zeta_i;
   mpzx_t factor, xplusa, pow, gcd;
   cm_timer_t clock, clock2;

   cm_timer_start (clock);

   if (f->deg <= 3)
      /* PARI implements the formula for degree 2, and, since version 2.15,
         also for degree 3. We may as well let it handle the case
         of degree 1. */
      cm_pari_oneroot (root, f, p);
   else {
      mpz_init (zeta);
      n = good_root_of_unity (zeta, p, f->deg);
      /* Fix a target degree of the factor for early abort to avoid more
         gcds when the factor is "small enough". The average degree of
         the gcd is f->deg / n; we stop at about twice that, with a bound
         guaranteed to be at most f->deg - 1 and at least 1 since
         2 <= n <= f->deg. */
      target = (2 * f->deg) / n - 1;
      if (debug)
         cm_file_printf ("    n = %i, target = %i\n", n, target);
      mpz_init (e);
      mpz_sub_ui (e, p, 1);
      mpz_divexact_ui (e, e, n);
      mpzx_init (xplusa, 1);
      mpz_set_ui (xplusa->coeff [1], 1);
      mpz_set_ui (xplusa->coeff [0], 0);
      mpzx_init (pow, f->deg - 1);
      mpzx_init (gcd, -1);
      mpz_init (zeta_i);
      mpzx_init (factor, -1);
      while (factor->deg == -1) {
         cm_timer_start (clock2);
         mpz_add_ui (xplusa->coeff [0], xplusa->coeff [0], 1);
         mpzx_pow_modmod (pow, xplusa, e, f, p);
         cm_timer_stop (clock2);
         if (debug)
            cm_file_printf ("    Time for power: %.1lf\n",
               cm_timer_get (clock2));
         mpz_set_ui (zeta_i, 1);
         if (pow->deg >= 1)
            for (i = 1;
               i <= n && (factor->deg == -1 || factor->deg > target);
               i++) {
               cm_timer_start (clock2);
               mpz_mul (zeta_i, zeta_i, zeta);
               mpz_mod (zeta_i, zeta_i, p); /* zeta^i */
               /* Shift the power and compute the gcd with f. */
               mpz_sub (pow->coeff [0], pow->coeff [0], zeta_i);
               mpz_mod (pow->coeff [0], pow->coeff [0], p);
               mpzx_gcd_mod (gcd, pow, f, p);
               /* Shift the power back. */
               mpz_add (pow->coeff [0], pow->coeff [0], zeta_i);
               mpz_mod (pow->coeff [0], pow->coeff [0], p);
               cm_timer_stop (clock2);
               if (debug)
                  cm_file_printf ("    Time for gcd, degree %i: %.1lf\n",
                     gcd->deg, cm_timer_get (clock2));
               if (gcd->deg >= 1) {
                  /* Consider the smaller one of gcd and f / gcd. Since gcd
                     usually has a low degree, this optimisation is of
                     interest only when f has low degree, so without much
                     impact overall. */
                  min = CM_MIN (gcd->deg, f->deg - gcd->deg);
                  if (factor->deg == -1 || min < factor->deg) {
                     if (min != gcd->deg)
                        mpzx_divexact_mod (gcd, f, gcd, p);
                     mpzx_clear (factor);
                     factor [0] = gcd [0];
                     mpzx_init (gcd, -1);
                  }
               }
            }
      }

      /* Recurse with the found factor. */
      mpzx_oneroot_split_mod_rec (root, factor, p, verbose, debug);
      mpz_clear (zeta);
      mpz_clear (e);
      mpzx_clear (xplusa);
      mpzx_clear (pow);
      mpzx_clear (gcd);
      mpz_clear (zeta_i);
      mpzx_clear (factor);
   }

   cm_timer_stop (clock);
}

/*****************************************************************************/

void mpzx_oneroot_split_mod (mpz_ptr root, mpzx_srcptr f, mpz_srcptr p,
   bool verbose, bool debug)
   /* Compute in root a root of the monic polynomial f over the prime field
      of characteristic p, assuming that f splits completely. */
{
   mpzx_t F;
   cm_timer_t clock;

   cm_timer_start (clock);
   if (verbose && f->deg > 1)
      cm_file_printf ("  Root finding in degree %i\n", f->deg);

   mpzx_init (F, f->deg);
   mpzx_mod (F, f, p);

   mpzx_oneroot_split_mod_rec (root, F, p, verbose, debug);

   mpzx_clear (F);

   cm_timer_stop (clock);
   if (verbose && f->deg > 1)
      cm_file_printf ("  Time for root: %.1f\n", cm_timer_get (clock));
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
   /* Return in prim the product of all primes p with a < p <= b.
      This function requires some memory to compute the prime list by
      a call to PARI. When calling the function, the PARI stack must be
      empty, so that it can be reset safely. */
{
   pari_sp av;
   GEN p;
   mpz_t *prod;
   int len, len2, i;

   av = avma;
   p = primes_interval_zv (a+1, b);
   len = glength (p);
   if (len == 0) {
      /* This probably never occurs... */
      len = 1;
      prod = (mpz_t *) malloc (1 * sizeof (mpz_t));
      mpz_init_set_ui (prod [0], 1);
   }
   else {
      prod = (mpz_t *) malloc (len * sizeof (mpz_t));
      for (i = 0; i < len; i++)
         mpz_init_set_si (prod [i], p [i+1]);
   }
   avma = av;
   /* The computations usually enlarge the PARI stack, but later steps
      do not need as much PARI memory, so it can be released. */
   parivstack_reset ();

   while (len > 1) {
      len2 = len / 2;
      for (i = 0; i < len2; i++)
         mpz_mul (prod [i], prod [2*i], prod [2*i+1]);
      if (len % 2 != 0) {
         mpz_set (prod [len2], prod [2*len2]);
         len2++;
      }
      for (i = len2; i < len; i++)
         mpz_clear (prod [i]);
      len = len2;
   }

   mpz_set (prim, prod [0]);
   mpz_clear (prod [0]);
   free (prod);
}

/*****************************************************************************/

bool cm_pari_cornacchia (mpz_ptr t, mpz_ptr v, mpz_srcptr p,
   mpz_srcptr root, const int_cl_t d)
{
   /* Compute t such that 4*p = t^2-v^2*d for some v, where p is an odd
      prime and d is an imaginary-quadratic discriminant such that d is a
      square modulo p and |d|<4*p.
      The return value indicates whether such a t exists; if not, the
      value of t is not changed during the algorithm. If yes and v is not
      NULL, it is changed.
      If root is not NULL, it is assumed to contain a pre-computed
      square root of d modulo p. */

   pari_sp av;
   GEN half, M, R;
   mpz_t rootloc;
   mpz_t r0, r1, ri, rim1, rip1;
      /* remainders 0, 1, i, i-1 and i+1 of the Euclidian algorithm */
   mpz_t tp, vp;
      /* candidates t' and v' for t and v */
   mpz_t l;
      /* stop Euclidian algorithm with ri > l >= rip1 */
   mpz_t Vi, Vip1, qi;
      /* Bezout coefficients in front of r1 at steps i and i+1,
         and quotient at step i */
   mpz_t tmp;
   bool ok;

   mpz_init (rootloc);
   mpz_init (r0);
   mpz_init (r1);
   mpz_init (rip1);
   mpz_init (tp);
   mpz_init (vp);
   mpz_init (tmp);

   /* Prepare a root of d modulo p. */
   if (root != NULL)
      mpz_set (rootloc, root);
   else
      cm_nt_mpz_tonelli_si (rootloc, d, p);

   /* Prepare the halfgcd. */
   if ((d - 1) % 8 == 0) {
      /* Solve p = (t/2)^2 + y^2 * (-d). */
      mpz_set (r0, p);
      mpz_set (r1, rootloc);
   }
   else if (d % 4 == 0) {
      /* Solve p = (t/2)^2 + y^2 * (-d/4); we need to divide root by 2. */
      mpz_set (r0, p);
      if (!mpz_divisible_2exp_p (rootloc, 1))
         mpz_sub (rootloc, p, rootloc);
      mpz_divexact_ui (r1, rootloc, 2);
   }
   else {
      /* d = 5 mod 8, the complicated case */
      mpz_init (l);
      mpz_init (ri);
      mpz_init (Vi);
      mpz_init (Vip1);

      mpz_mul_2exp (r0, p, 1);
      /* Make root odd, then it is a root of d modulo 4*p. */
      if (mpz_divisible_2exp_p (rootloc, 1))
         mpz_sub (r1, p, rootloc);
      else
         mpz_set (r1, rootloc);
      mpz_mul_2exp (l, p, 2);
      mpz_sqrt (l, l);
   }

   /* Delegate the halfgcd to pari. */
   av = avma;
   half = ghalfgcd (mpz_get_Z (r0), mpz_get_Z (r1));
   R = gel (half, 2);
   Z_get_mpz (rip1, gel (R, 2));
   if ((d - 5) % 8 == 0) {
      Z_get_mpz (ri, gel (R, 1));
      M = gel (half, 1);
      Z_get_mpz (Vi, gcoeff (M, 1, 2));
      Z_get_mpz (Vip1, gcoeff (M, 2, 2));
   }
   avma = av;

   /* Determine the candidate for t. */
   if ((d - 5) % 8 != 0)
      mpz_mul_2exp (tp, rip1, 1);
   else {
      if (mpz_cmp (ri, l) > 0)
         mpz_set (tp, rip1);
      else {
         /* Compute the previous remainder r_{i-1}. */
         mpz_init (rim1);
         mpz_init (qi);
         mpz_tdiv_q (qi, Vip1, Vi);
         mpz_abs (qi, qi);
         mpz_mul (rim1, ri, qi);
         mpz_add (rim1, rim1, rip1);

         if (mpz_cmp (rim1, l) > 0)
            mpz_set (tp, ri);
         else
            /* Now we have sqrt (4*p) > r_{i-1} > r_i > sqrt (2*p)
               and r_{i+1} < sqrt (2*p), since by the properties of the
               halfgcd r_{i+1} is the first remainder below this bound.
               This implies
               r_{i-2} = q_{i-1} * r_{i-1}    + r_i
                       >     1   * sqrt (2*p) + sqrt (2*p)
                       > sqrt (4*p),
               so there is no need to go further up. */
            mpz_set (tp, rim1);

         mpz_clear (rim1);
         mpz_clear (qi);
      }

      mpz_clear (l);
      mpz_clear (ri);
      mpz_clear (Vi);
      mpz_clear (Vip1);
   }

   /* Check whether v exists. */
   mpz_mul_2exp (vp, p, 2);
   mpz_pow_ui (tmp, tp, 2);
   mpz_sub (vp, vp, tmp);
   if (!mpz_divisible_ui_p (vp, -d))
      ok = false;
   else {
      mpz_divexact_ui (vp, vp, -d);
      if (!mpz_perfect_square_p (vp))
            ok = false;
      else
            ok = true;
   }
   if (ok) {
      mpz_set (t, tp);
      if (v != NULL)
         mpz_root (v, vp, 2);
   }

   mpz_clear (rootloc);
   mpz_clear (r0);
   mpz_clear (r1);
   mpz_clear (rip1);
   mpz_clear (tp);
   mpz_clear (vp);
   mpz_clear (tmp);

   return ok;
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
