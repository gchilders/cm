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
#include <pari/pari.h>
#include "cm_class-impl.h"

static GEN mpz_get_Z (mpz_t z);
static void Z_get_mpz (mpz_t z, GEN x);
static GEN mpzx_get_FpX (mpz_t *f, int deg, mpz_t p);
static int ZX_get_mpzx (mpz_t *res, GEN f);
static void cm_pari_onefactor (mpz_t *res, mpz_t *f, int deg, int deg_factor,
   mpz_t p, bool verbose);

/*****************************************************************************/
/*                                                                           */
/* Conversion functions between PARI, GMP and internal types.                */
/*                                                                           */
/*****************************************************************************/

static GEN mpz_get_Z (mpz_t z)
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

static void Z_get_mpz (mpz_t z, GEN x)
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

static GEN mpzx_get_FpX (mpz_t *f, int deg, mpz_t p)
   /* returns a GEN of type t_POL over t_INT corresponding to the polynomial */
   /* f of degree deg, reduced modulo p                                      */
   /* The zero polynomial is characterised by a degree equal to -1.          */

{
   int i;
   GEN res;
   mpz_t tmp;

   mpz_init (tmp);

   res = cgetg (deg + 3, t_POL);
   setvarn (res, 0);
   if (deg == -1)
      setsigne (res, 0);
   else {
      setsigne (res, 1);
      for (i = 0; i <= deg; i++) {
         mpz_mod (tmp, f [i], p);
         gel (res, i+2) = mpz_get_Z (tmp);
      }
   }

   mpz_clear (tmp);

   return res;
}

/*****************************************************************************/

static int ZX_get_mpzx (mpz_t *res, GEN f)
   /* f is a GEN of type t_POL over t_INT                                    */
   /* The function returns the degree of f and, via res, the array of its    */
   /* coefficients. res must contain sufficient space to hold all these, and */
   /* its entries must already be initialised with mpz_init.                 */

{
   int i;
   int deg = degree (f);

   for (i = 0; i <= deg; i++)
      Z_get_mpz (res [i], gel (f, i+2));

   return deg;
}

/*****************************************************************************/
/*                                                                           */
/* Functions for finding factors of polynomials.                             */
/*                                                                           */
/*****************************************************************************/

static GEN good_root_of_unity (int *n, const GEN p, const int deg,
   const int deg_factor)
   /* returns a root of unity in F_p that is suitable for finding a factor   */
   /* of degree deg_factor of a polynomial of degree deg; the order is       */
   /* returned in n                                                          */
   /* A good choice seems to be n close to deg/deg_factor; we choose n       */
   /* twice as big and decrement until it divides p-1.                       */

{
   GEN pm, factn, power, base, zeta;

   pari_sp lbot, ltop = avma;

   pm = subis (p, 1ul);
   for (*n = deg / 2 / deg_factor + 1; !dvdiu (pm, *n); (*n)--);

   factn = Z_factor (stoi (*n));
   power = diviuexact (pm, *n);
   base = gen_1;
   do {
      base = addis (base, 1l);
      zeta = Fp_pow (base, power, p);
      lbot = avma; /* memorise the only interesting object */
   }
   while (!equaliu (Fp_order (zeta, factn, p), *n));

   avma = lbot;
   zeta = gerepileuptoint (ltop, zeta);
   return zeta;
}

/*****************************************************************************/

static void cm_pari_onefactor (mpz_t *res, mpz_t *f, int deg, int deg_factor,
   mpz_t p, bool verbose)
   /* Assuming that deg_factor is the minimal degree of a factor of f, finds */
   /* such a factor and returns it in res. The coefficients of res need to   */
   /* already be initialised.                                                */

{
   GEN pp;
   int n, deg_f, target, i;
   GEN prim, expo, fact, minfactor, xplusa, zeta, xpow, tmp;
   cm_timer clock, clock2;
   pari_sp av;

   cm_timer_start (clock);
   if (verbose)
      printf ("--- Factor finding, degree %i out of %i\n", deg_factor, deg);

   pari_init (500000, 0);
   paristack_setsize (500000, 500000000);
   av = avma;

   pp = mpz_get_Z (p);
   fact = mpzx_get_FpX (f, deg, p);
   minfactor = fact; /* factor of minimal degree found so far */

   xplusa = pol_x (varn (fact));
   while (degpol (minfactor) != deg_factor) {
      /* split minfactor by computing its gcd with (X+a)^exp-zeta, where
         zeta varies over the n-th roots of unity in F_p for a suitable n */
      fact = FpX_normalize (minfactor, pp);
      deg_f = degpol (fact);
      prim = good_root_of_unity (&n, pp, deg_f, deg_factor);
      expo = diviuexact (subis (powiu (pp, (unsigned long int) deg_factor), 1),
                         (unsigned long int) n);
      /* update X+a, avoid a=0 */
      gel (xplusa, 2) = addis (gel (xplusa, 2), 1);
      cm_timer_start (clock2);
      xpow = FpXQ_pow (xplusa, expo, fact, pp);
      cm_timer_stop (clock2);
      if (verbose)
         printf ("- Time for pow with n=%i in degree %i: %.1f\n",
                 n, deg_f, cm_timer_get (clock2));
      if (degpol (xpow) > 0 && degpol (xpow) < deg_f) {
         /* Fix a target for early abort, depending on n and deg_f, to
            avoid computing too many gcds; we stop as soon as a degree of
            at most target is reached. It is necessary to have
            target <= deg_f - 1. */
         target = (2 * deg_f) / n - 1;
         for (i = 0, zeta = gen_1;
            i < n
               && degree (minfactor) > deg_factor
               /* stop early to avoid too many gcds */
               && degree (minfactor) > target;
            i++, zeta = Fp_mul (zeta, prim, pp)) {
            cm_timer_start (clock2);
            tmp = FpX_gcd (FpX_Fp_sub (xpow, zeta, pp), fact, pp);
            cm_timer_stop (clock2);
            if (verbose)
               printf ("- Time for gcd of degrees %li, %li -> %li: %.1f\n",
                  degpol (fact), degpol (xpow), degpol (tmp),
                  cm_timer_get (clock2));
            if (degpol (tmp) > 0 && degpol (tmp) < degpol (fact)) {
               fact = FpX_div (fact, tmp, pp);
               if (degpol (tmp) < degpol (minfactor))
                  minfactor = tmp;
            }
         }
      }
    }

   ZX_get_mpzx (res, FpX_normalize (minfactor, pp));

   avma = av;
   pari_close ();

   cm_timer_stop (clock);
   if (verbose)
      printf ("-- Time for factor: %.1f\n", cm_timer_get (clock));
}

/*****************************************************************************/

#if 0
void cm_pari_oneroot_alt (mpz_t root, mpz_t *f, int deg, mpz_t p, bool verbose)
   /* finds a root of the polynomial f of degree deg over the prime field of */
   /* characteristic p, assuming that such a root exists, and returns it in  */
   /* the variable of the same name                                          */

{
   GEN fp, pp, rootp;
   cm_timer clock;

   cm_timer_start (clock);
   if (verbose)
      printf ("--- Root finding in degree %i\n", deg);

   pari_init (2000 * deg * mpz_sizeinbase (p, 2) / 8 + 1000000, 0);

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, deg, p);

   rootp = FpX_oneroot (fp, pp);
   Z_get_mpz (root, rootp);

   pari_close ();

   cm_timer_stop (clock);
   if (verbose)
      printf ("-- Time for root: %.1f\n", cm_timer_get (clock));
}
#endif

/*****************************************************************************/

void cm_pari_oneroot (mpz_t root, mpz_t *f, int deg, mpz_t p, bool verbose)
   /* finds a root of the polynomial f of degree deg over the prime field of */
   /* characteristic p, assuming that such a root exists, and returns it in  */
   /* the variable of the same name                                          */

{
   mpz_t fact [2];

   mpz_init (fact [0]);
   mpz_init (fact [1]);
   cm_pari_onefactor (fact, f, deg, 1, p, verbose);
   if (mpz_sgn (fact [0]) != 0)
      mpz_sub (root, p, fact [0]);
   else
      mpz_set_ui (root, 0ul);
   mpz_clear (fact [0]);
   mpz_clear (fact [1]);
}

/*****************************************************************************/

mpz_t* cm_pari_find_roots (mpz_t *f, int deg, mpz_t p, int *no)
   /* computes all the roots of the polynomial f of degree deg in the prime  */
   /* field of characteristic p without multiplicities                       */
   /* no is the number of found roots                                        */

{
   mpz_t *res;
   GEN fp, pp, rootsp;
   int i;

   pari_init (2000 * deg * mpz_sizeinbase (p, 2) / 8 + 1000000, 0);

   pp = mpz_get_Z (p);
   fp = mpzx_get_FpX (f, deg, p);
   rootsp = FpX_roots (fp, pp);
   *no = lg (rootsp) - 1;
   res = (mpz_t*) malloc ((*no) * sizeof (mpz_t));
   for (i = 0; i < *no; i++) {
      mpz_init (res [i]);
      Z_get_mpz (res [i], gel (rootsp, i+1));
   }

   pari_close ();

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

   pari_init (500000, 0);
   paristack_setsize (500000, 500000000);
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
   pari_close ();

   return length;
}


/*****************************************************************************/
