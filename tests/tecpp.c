/*

tecpp.c - test of ECPP

Copyright (C) 2021, 2022 Andreas Enge

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

static void test_ecpp (mpz_srcptr n)
   /* Try to compute an ECPP certificate for n. */
{
   bool res;

   res =    cm_ecpp (n, CM_MODPOLDIR,
      NULL /* filename */,
      NULL /* tmpdir */,
      false /* print */,
      true /* trust */,
      true /* check */,
      false /* verbose */,
      false /* debug */);

   if (!res) {
      printf ("ECPP certificate verification failed for\n");
      mpz_out_str (stdout, 10, n);
      printf ("\n");
      exit (1);
   }
}

/*****************************************************************************/

int main (void)
{
   mpz_t n;

   mpz_init (n);
   mpz_set_ui (n, 10);
   mpz_pow_ui (n, n, 200);
   mpz_nextprime (n, n);

   cm_pari_init ();
   test_ecpp (n);
   cm_pari_clear ();

   mpz_clear (n);

   return 0;
}
