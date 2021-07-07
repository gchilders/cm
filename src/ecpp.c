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

   /* The following are annoying numbers for which it is difficult to
      find suitable discriminants. */
//   mpz_set_str (N, "6921429791794771459742425597804321859721736769688598662256548209942641251770313333951880454803733175805128147717974272402897979725899194689808268585820286779821073226454686860838504940028382260983619367446243748011228298985156480351736074698424386132694515069696917811709986474707015305812933093675728793595641150023900317151001253198621596601096017878695667648918896345968275829768823240980678309332041538791440704304588270662770485615405711160624896327085530458306558095928233591221254691215172835594088919257205266254956364078291052319835486199828670449883501465041676641589250575498569088624307472553053606299442533281428150007691854115634707257865200910074176320902009918653319891314929141045724434416043992175569308334667825431094630024816975158767932366431318707237808241231387597850283688855754437837939066389526337543833786229458437144045180777406933260774887340771999870577010956023", 10);

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
