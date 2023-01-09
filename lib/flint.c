/*

flint.c - functions using flint; for finding roots of polynomials

Copyright (C) 2022 Andreas Enge

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

#endif

