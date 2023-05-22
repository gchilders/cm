/*

flint.c - functions using flint

Copyright (C) 2022, 2023 Andreas Enge

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

#ifdef HAVE_FLINT

/*****************************************************************************/
/*                                                                           */
/* Conversion functions between FLINT and internal types.                    */
/*                                                                           */
/*****************************************************************************/

void mpzx_set_fmpz_mod_poly (mpzx_ptr f, fmpz_mod_poly_t ff,
   const fmpz_mod_ctx_t ctx)
{
   int deg, i;

   deg = fmpz_mod_poly_degree (ff, ctx);
   mpzx_set_deg (f, deg);
   for (i = 0; i <= deg; i++)
      fmpz_mod_poly_get_coeff_mpz (f->coeff [i], ff, i, ctx);
}

/*****************************************************************************/

void fmpz_mod_poly_set_mpzx (fmpz_mod_poly_t ff, mpzx_srcptr f,
   const fmpz_mod_ctx_t ctx)
{
   int deg, i;

   deg = f->deg;
   fmpz_mod_poly_realloc (ff, deg + 1, ctx);
   for (i = 0; i <= deg; i++)
      fmpz_mod_poly_set_coeff_mpz (ff, i, f->coeff [i], ctx);
}

/*****************************************************************************/
/*                                                                           */
/* Functions for mpzx modulo p relying on FLINT.                             */
/*                                                                           */
/*****************************************************************************/

void cm_flint_mpzx_xplusa_pow_modmod (mpzx_ptr g, unsigned long int a,
   mpz_srcptr e, mpzx_srcptr m, mpz_srcptr p)
   /* Compute g = (X+a)^e modulo m and p. */
{
   fmpz_t pp, ep, ap;
   fmpz_mod_ctx_t ctx;
   fmpz_mod_poly_t mp, gp, minv;

   fmpz_init (pp);
   fmpz_set_mpz (pp, p);
   fmpz_init (ap);
   fmpz_init (ep);
   fmpz_mod_ctx_init (ctx, pp);
   fmpz_mod_poly_init (mp, ctx);
   fmpz_mod_poly_init (gp, ctx);
   fmpz_mod_poly_init (minv, ctx);

   fmpz_set_mpz (ep, e);
   fmpz_set_ui (ap, a);
   fmpz_mod_poly_set_mpzx (mp, m, ctx);
   fmpz_mod_poly_reverse (minv, mp, mp->length, ctx);
   fmpz_mod_poly_inv_series (minv, minv, mp->length, ctx);

   fmpz_mod_poly_powmod_linear_fmpz_preinv (gp, ap, ep, mp, minv, ctx);

   mpzx_set_fmpz_mod_poly (g, gp, ctx);

   fmpz_clear (pp);
   fmpz_clear (ep);
   fmpz_clear (ap);
   fmpz_mod_poly_clear (mp, ctx);
   fmpz_mod_poly_clear (gp, ctx);
   fmpz_mod_poly_clear (minv, ctx);
   fmpz_mod_ctx_clear (ctx);
}

/*****************************************************************************/

void cm_flint_mpzx_gcd_mod (mpzx_ptr h, mpzx_srcptr f, mpzx_srcptr g,
   mpz_srcptr p)
   /* Compute h = gcd (f, g) modulo p without imposing a normalisation of the
      leading coefficient of h. */
{
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
}

#endif

