/*

pari.c - functions for factoring polynomials using pari

Copyright (C) 2010, 2015 Andreas Enge

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
#include <pari/pari.h>
#include "cm_class-impl.h"

static GEN mpz_get_Z (mpz_t z);
static void Z_get_mpz (mpz_t z, GEN x);
static GEN mpzx_get_FpX (mpz_t *f, int deg, mpz_t p);
static int ZX_get_mpzx (mpz_t *res, GEN f);

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

static GEN good_root_of_unity (int *n, const GEN p, const int deg,
   const int deg_factor, bool verbose)
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
   if (verbose)
      printf ("n %i\n", *n);

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

void cm_pari_onefactor (mpz_t *res, mpz_t *f, int deg, int deg_factor,
   mpz_t p, bool verbose)
   /* Assuming that deg_factor is the minimal degree of a factor of f, finds */
   /* such a factor and returns it in res. The coefficients of res need to   */
   /* already be initialised.                                                */

{
   GEN pp;
   int n, deg_f, i;
   GEN prim, expo, fact, minfactor, xplusa, zeta, xpow, tmp;
   cm_timer clock, clock2;

   cm_timer_start (clock);
   if (verbose)
      printf ("--- Factor finding, degree %i out of %i\n", deg_factor, deg);

   pari_init (2000 * deg * mpz_sizeinbase (p, 2) / 8 + 1000000, 0);

   pp = mpz_get_Z (p);
   fact = mpzx_get_FpX (f, deg, p);
   minfactor = fact; /* factor of minimal degree found so far */

   prim = good_root_of_unity (&n,pp, deg, deg_factor, verbose);
   expo = diviuexact (subis (powiu (pp, (unsigned long int) deg_factor), 1),
            (unsigned long int) n);

   xplusa = pol_x (varn (fact));
   zeta = gen_1;
   while (degpol (minfactor) != deg_factor) {
      /* split minfactor by computing its gcd with (X+a)^exp-zeta, where    */
      /* zeta varies over the roots of unity in F_p                         */
      fact = FpX_normalize (minfactor, pp);
      deg_f = degpol (fact);
      /* update X+a, avoid a=0 */
      gel (xplusa, 2) = addis (gel (xplusa, 2), 1);
      cm_timer_start (clock2);
      xpow = FpXQ_pow (xplusa, expo, fact, pp);
      cm_timer_stop (clock2);
      if (verbose)
         printf ("- Time for pow: %.1f in degree %i\n", cm_timer_get (clock2),
                 deg_f);
      for (i = 0; i < n; i++) {
         cm_timer_start (clock2);
         tmp = FpX_gcd (FpX_Fp_sub (xpow, zeta, pp), fact, pp);
         cm_timer_stop (clock2);
         if (verbose)
            printf ("- Time for gcd of degrees %li, %li -> %li: %.1f\n",
               degpol (fact), degpol (xpow), degpol (tmp),
               cm_timer_get (clock2));
         if (degpol (tmp) > 0 && degpol (tmp) < degpol (fact)) {
            fact = FpX_div (fact, tmp, pp);
            if (degpol (tmp) < degpol (minfactor)) {
               minfactor = tmp;
               if (degree (minfactor) == deg_factor ||
                   degree (minfactor) <= deg_f / (n / 2) + 1)
                  /* stop early to avoid too many gcds */
                  break;
            }
         }
         zeta = Fp_mul (zeta, prim, pp);
      }
    }

   ZX_get_mpzx (res, FpX_normalize (minfactor, pp));

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
