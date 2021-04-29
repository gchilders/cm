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
/*                                                                           */
/* Functions for finding factors of polynomials.                             */
/*                                                                           */
/*****************************************************************************/

void cm_pari_oneroot (mpz_t root, mpz_t *f, int deg, mpz_t p, bool verbose)
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

   rootp = FpX_oneroot_split (fp, pp);
   Z_get_mpz (root, rootp);

   pari_close ();

   cm_timer_stop (clock);
   if (verbose)
      printf ("-- Time for root: %.1f\n", cm_timer_get (clock));
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
