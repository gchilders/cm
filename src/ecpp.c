/*

ecpp.c - executable for elliptic curve primality proofs

Copyright (C) 2021 Andreas Enge

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

#include "params.h"

int main (void)
{
   mpz_t N;

   cm_pari_init ();
   mpz_init (N);

   mpz_set_ui (N, 10);
   mpz_pow_ui (N, N, 1000);
   mpz_nextprime (N, N);

   cm_ecpp (N, CM_MODPOLDIR,
      false /* pari */,
      true /* tower */,
      false /* print */,
      true /* verbose */,
      true /* debug */);

   cm_pari_clear ();
   mpz_clear (N);

   return 0;
}
