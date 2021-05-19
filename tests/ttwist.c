/*

tcm.c - tests with twisted curves

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

#include "cm.h"

/*****************************************************************************/

static void test_curve (int_cl_t d, mpz_srcptr p, mpz_srcptr l,
   mpz_srcptr co)
   /* Try to find a curve of discriminant d over Fp with a point of
      order l and cofactor co, where co is coprime to l (and l is not
      necessarily prime). */
{
   cm_param_t param;
   cm_class_t c;
   mpz_t a, b, x, y;

   mpz_init (a);
   mpz_init (b);
   mpz_init (x);
   mpz_init (y);
   
   printf ("p="); mpz_out_str (stdout, 10, p); printf ("\n");
   printf ("l="); mpz_out_str (stdout, 10, l); printf ("\n");
   printf ("c="); mpz_out_str (stdout, 10, co); printf ("\n");
   cm_param_init (param, d, CM_INVARIANT_J, false);
   cm_class_init (c, param, false);
   cm_class_compute (c, param, false, true, false);
   cm_curve_and_point (a, b, x, y, param, c, p, l, co, CM_MODPOLDIR,
      false);
      /* CM_MODPOLDIR is a preprocessor variable defined by the -D
         parameter of gcc */
   cm_class_clear (c);
   printf ("ok\n");

   mpz_clear (a);
   mpz_clear (b);
   mpz_clear (x);
   mpz_clear (y);
}

/*****************************************************************************/

int main (void)
{
   mpz_t p, l, c;

   cm_pari_init ();
   mpz_init (p);
   mpz_init (l);
   mpz_init (c);

   /* Test whether the code correctly constructs all four twists for D=-4.
      With the following p, 2 generates Fp* and in particular Fp* modulo
      fourth powers, and the four twists are given by a=1, 2, 4, 8 and
      b = 0. Their cardinalities are easily computed by PARI/GP, with the
      results as follows. */
   mpz_set_str (p, "115792089237316195398462578067141184803371344843544913790019434091196253893701", 10);
   mpz_set_str (l, "1157920892373161953984625780671411848033713448435449137900194340911962538929", 10);
   mpz_set_ui (c, 10);
   test_curve (-4, p, l, c);

   mpz_set_str (l, "57896044618658097699231289033570592402025954788693395358436198164882476055081", 10);
   mpz_set_ui (c, 2);
   test_curve (-4, p, l, c);

   mpz_set_str (l, "14474011154664524424807822258392648100421418105443114223752429261399531736813", 10);
   mpz_set_ui (c, 8);
   test_curve (-4, p, l, c);

   mpz_set_str (l, "57896044618658097699231289033570592401345390054851518431583235926313777838621", 10);
   mpz_set_ui (c, 2);
   test_curve (-4, p, l, c);

   /* Test whether the code correctly constructs all six twists for
      D=-3. With the following p, 17 is the smallest generator of Fp* and
      also of Fp* modulo sixth powers. */
   mpz_set_str (p, "1606938044258984566551191268523188429582205066354293683213377", 10);
   mpz_set_str (l, "44637167896082904626421979681234891116178702597333071464893", 10);
   mpz_set_ui (c, 36);
   test_curve (-3, p, l, c);

   mpz_set_str (l, "1606938044258984566551191268521920778981976839204596793695879", 10);
   mpz_set_ui (c, 1);
   test_curve (-3, p, l, c);

   mpz_set_str (l, "535646014752994855517063756173551042793916204018299968057703", 10);
   mpz_set_ui (c, 3);
   test_curve (-3, p, l, c);

   mpz_set_str (l, "200867255532373070818898908565240097372747104900574599211326", 10);
   mpz_set_ui (c, 8);
   test_curve (-3, p, l, c);

   cm_pari_clear ();
   mpz_clear (p);
   mpz_clear (l);
   mpz_clear (c);

   return 0;
}
