/*

jmodp.c - code for obtaining a j-invariant from a class polynomial

Copyright (C) 2009, 2010, 2021 Andreas Enge

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

static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose);
static mpz_t* get_j_mod_P_from_modular (int *no, const char* modpoldir,
   char type, int level, mpz_t root, mpz_t P);
static mpz_t* simpleeta_cm_get_j_mod_P (cm_class_t c, mpz_t root, mpz_t P,
   int *no);

/*****************************************************************************/

static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose)
   /* returns a root of the minimal polynomial modulo P                      */
   /* root is changed                                                        */

{
   /* The local variables are needed only for the complex case. */
   cm_timer clock;
   mpzx_t minpoly_P;
   int i;
   mpz_t omega, tmp;
   /* the second element of the integral basis of the maximal order */
   /* modulo P                                                      */

   if (c.field == CM_FIELD_REAL)
      cm_pari_oneroot (root, c.minpoly, P, verbose);
   else {
      mpz_init (omega);
      mpz_init (tmp);

      mpzx_init (minpoly_P, c.minpoly->deg);

      /* compute the second element of the integral basis modulo P */
      cm_timer_start (clock);
      cm_nt_mpz_tonelli (omega, c.dfund, P);
      cm_timer_stop (clock);
      if (verbose)
         printf ("--- Time for square root: %.1f\n", cm_timer_get (clock));
      if (c.dfund % 4 != 0)
         mpz_add_ui (omega, omega, 1ul);
      mpz_set_ui (tmp, 2ul);
      mpz_invert (tmp, tmp, P);
      mpz_mul (omega, omega, tmp);
      mpz_mod (omega, omega, P);

      /* create the minimal polynomial modulo P */
      for (i = 0; i < c.minpoly->deg; i++) {
         mpz_mod (minpoly_P->coeff [i], c.minpoly_complex->coeff [i], P);
         mpz_mul (tmp, minpoly_P->coeff [i], omega);
         mpz_mod (minpoly_P->coeff [i], tmp, P);
         mpz_add (tmp, minpoly_P->coeff [i], c.minpoly->coeff [i]);
         mpz_mod (minpoly_P->coeff [i], tmp, P);
      }
      mpz_set_ui (minpoly_P->coeff [c.minpoly->deg], 1ul);

      cm_pari_oneroot (root, minpoly_P, P, verbose);

      mpz_clear (omega);
      mpz_clear (tmp);
      mpzx_clear (minpoly_P);
   }

   if (verbose) {
      printf ("Root: ");
      mpz_out_str (stdout, 10, root);
      printf ("\n");
   }
}

/*****************************************************************************/

static mpz_t* get_j_mod_P_from_modular (int *no, const char* modpoldir,
   char type, int level, mpz_t root, mpz_t P)
   /* computes the possible j-values as roots of the modular polynomial      */
   /* specified by type and level                                            */

{
   mpz_t  *j;
   mpz_t* poly_j;
      /* a polynomial one of whose roots is j mod P */
   int    n, i;

   poly_j = cm_modpol_read_specialised_mod (&n, level, type, P, root,
      modpoldir);
   j = cm_pari_find_roots (poly_j, n, P, no);
   for (i = 0; i <= n; i++)
      mpz_clear (poly_j [i]);
   free (poly_j);

   return j;
}

/*****************************************************************************/

static mpz_t* simpleeta_cm_get_j_mod_P (cm_class_t c, mpz_t root, mpz_t P,
   int *no)

{
   mpz_t* j = (mpz_t*) malloc (sizeof (mpz_t));
   mpz_t f3, tmp;

   mpz_init (j [0]);
   mpz_init (f3);
   mpz_init (tmp);

   mpz_powm_ui (root, root, (unsigned long int) (c.s / c.e), P);

   if (c.p [0] == 3) {
      mpz_add_ui (f3, root, 3ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_invert (tmp, root, P);
      mpz_mul_ui (tmp, tmp, 27ul);
      mpz_add_ui (tmp, tmp, 1ul);
      mpz_mul (j [0], f3, tmp);
      mpz_mod (j [0], j [0], P);
   }
   else if (c.p [0] == 5) {
      mpz_add_ui (tmp, root, 10ul);
      mpz_mul (f3, root, tmp);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 5ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_invert (tmp, root, P);
      mpz_mul (j [0], f3, tmp);
      mpz_mod (j [0], j [0], P);
   }
   else if (c.p [0] == 7) {
      mpz_add_ui (tmp, root, 5ul);
      mpz_mul (f3, root, tmp);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 1ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_add_ui (j [0], root, 13ul);
      mpz_mul (j [0], j [0], root);
      mpz_mod (j [0], j [0], P);
      mpz_add_ui (j [0], j [0], 49ul);
      mpz_mul (j [0], j [0], f3);
      mpz_mod (j [0], j [0], P);
      mpz_invert (tmp, root, P);
      mpz_mul (j [0], j [0], tmp);
      mpz_mod (j [0], j [0], P);
   }
   else if (c.p [0] == 13) {
      mpz_add_ui (f3, root, 7ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 20ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 19ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 1ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_add_ui (j [0], root, 5ul);
      mpz_mul (j [0], j [0], root);
      mpz_mod (j [0], j [0], P);
      mpz_add_ui (j [0], j [0], 13ul);
      mpz_mul (j [0], j [0], f3);
      mpz_mod (j [0], j [0], P);
      mpz_invert (tmp, root, P);
      mpz_mul (j [0], j [0], tmp);
      mpz_mod (j [0], j [0], P);
   }
   else if (c.p [0] == 4) {
      mpz_add_ui (f3, root, 16ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_invert (tmp, f3, P);
      mpz_add_ui (f3, f3, 16ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_mul (j [0], tmp, f3);
      mpz_mod (j [0], j [0], P);
   }
   else if (c.p [0] == 9) {
      mpz_add_ui (tmp, root, 9ul);
      mpz_mul (tmp, tmp, root);
      mpz_mod (tmp, tmp, P);
      mpz_add_ui (tmp, tmp, 27ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_invert (j [0], tmp, P);
      mpz_add_ui (f3, tmp, 3ul);
      mpz_add_ui (tmp, root, 3ul);
      mpz_mul (f3, f3, tmp);
      mpz_mod (f3, f3, P);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_mul (j [0], j [0], f3);
      mpz_mod (j [0], j [0], P);
   }
   else if (c.p [0] == 25) {
      mpz_add_ui (f3, root, 10ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 55ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 200ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 525ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 1010ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 1425ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 1400ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 875ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 250ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 5ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_add_ui (tmp, root, 5ul);
      mpz_mul (tmp, tmp, root);
      mpz_mod (tmp, tmp, P);
      mpz_add_ui (tmp, tmp, 15ul);
      mpz_mul (tmp, tmp, root);
      mpz_mod (tmp, tmp, P);
      mpz_add_ui (tmp, tmp, 25ul);
      mpz_mul (tmp, tmp, root);
      mpz_mod (tmp, tmp, P);
      mpz_add_ui (tmp, tmp, 25ul);
      mpz_mul (tmp, tmp, root);
      mpz_mod (tmp, tmp, P);
      mpz_invert (j [0], tmp, P);
      mpz_mul (j [0], j [0], f3);
      mpz_mod (j [0], j [0], P);
   }

   *no = 1;

   mpz_clear (f3);
   mpz_clear (tmp);

   return j;
}

/*****************************************************************************/

mpz_t* cm_class_get_j_mod_P (int_cl_t d, char inv, mpz_t P, int *no,
   const char* modpoldir, bool readwrite, bool verbose)
   /* If readwrite is true, tries to read the polynomial from a file, and    */
   /* computes it and writes it to the file if it is not yet present.        */

{
   cm_class_t c;
   mpz_t *j;
   mpz_t root, d_mpz, tmp, tmp2, f24;
   cm_timer clock;

   cm_class_init (&c, d, inv, verbose);
   if (!readwrite || !cm_class_read (c))
      cm_class_compute_minpoly (c, false, readwrite, false, verbose);

   cm_timer_start (clock);
   mpz_init (root);
   get_root_mod_P (c, root, P, verbose);
   switch (inv)
   {
      case CM_INVARIANT_J:
         j = (mpz_t*) malloc (sizeof (mpz_t));
         mpz_init_set (j [0], root);
         *no = 1;
         break;
      case CM_INVARIANT_GAMMA2:
         j = (mpz_t*) malloc (sizeof (mpz_t));
         mpz_init_set (j [0], root);
         mpz_powm_ui (j [0], j [0], 3, P);
         *no = 1;
         break;
      case CM_INVARIANT_GAMMA3:
         j = (mpz_t*) malloc (sizeof (mpz_t));
         mpz_init (j [0]);

         mpz_init_set_si (d_mpz, d);
         mpz_init (tmp);
         mpz_init (tmp2);

         mpz_powm_ui (tmp, root, 2, P);
         mpz_invert (tmp2, d_mpz, P);
         mpz_mul (root, tmp, tmp2);
         mpz_add_ui (root, root, 1728);
         mpz_mod (j [0], root, P);

         *no = 1;

         mpz_clear (d_mpz);
         mpz_clear (tmp);
         mpz_clear (tmp2);
         break;
      case CM_INVARIANT_WEBER:
         j = (mpz_t*) malloc (sizeof (mpz_t));

         mpz_init (j [0]);
         mpz_init (f24);
         mpz_init (tmp);

         if (c.d % 3 == 0)
            mpz_powm_ui (f24, root, 2ul, P);
         else
            mpz_powm_ui (f24, root, 6ul, P);

         if (c.p [0] == 1) {
            mpz_mul_2exp (tmp, f24, 3ul);
            mpz_mod (tmp, tmp, P);
            mpz_powm_ui (f24, tmp, 2ul, P);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p [0] == 3) {
            mpz_powm_ui (tmp, f24, 4ul, P);
            mpz_set (f24, tmp);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p [0] == 5) {
            mpz_mul_2exp (tmp, f24, 6ul);
            mpz_set (f24, tmp);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p [0] == 7) {
            mpz_mul_2exp (tmp, f24, 3ul);
            mpz_mod (tmp, tmp, P);
            mpz_powm_ui (f24, tmp, 4ul, P);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p [0] == 2 || c.p [0] == 6) {
            mpz_mul_2exp (tmp, f24, 3ul);
            mpz_mod (tmp, tmp, P);
            mpz_powm_ui (f24, tmp, 2ul, P);
            mpz_add_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else {
            /* c.p [0] == 4 */
            mpz_mul_2exp (tmp, f24, 9ul);
            mpz_set (f24, tmp);
            mpz_add_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }

         *no = 1;

         mpz_clear (f24);
         mpz_clear (tmp);
         break;
      case CM_INVARIANT_DOUBLEETA:
#if 0
      case CM_INVARIANT_RAMIFIED:
#endif
         if (c.s != c.e)
            mpz_powm_ui (root, root, (unsigned long int) (c.s / c.e), P);
         j = get_j_mod_P_from_modular (no, modpoldir, CM_MODPOL_DOUBLEETA,
            c.p [0] * c.p [1], root, P);
         break;
      case CM_INVARIANT_MULTIETA:
      {
         int N = 1, i;
         for (i = 0; c.p [i] != 0; i++)
            N *= c.p [i];
         if (c.s != c.e)
            mpz_powm_ui (root, root, (unsigned long int) (c.s / c.e), P);
         j = get_j_mod_P_from_modular (no, modpoldir,
            CM_MODPOL_MULTIETA, N, root, P);
         break;
      }
      case CM_INVARIANT_ATKIN:
         j = get_j_mod_P_from_modular (no, modpoldir, CM_MODPOL_ATKIN,
            c.p [0], root, P);
         break;
      case CM_INVARIANT_SIMPLEETA:
         j = simpleeta_cm_get_j_mod_P (c, root, P, no);
         break;
      default: /* should not occur */
         printf ("class_cm_get_j_mod_P called for unknown class ");
         printf ("invariant\n");
         exit (1);
   }
   mpz_clear (root);
   cm_timer_stop (clock);
   if (verbose) {
      int i;

      printf ("%i candidate", *no);
      if (*no > 1)
         printf ("s");
      printf (" for j:");
      for (i = 0; i < *no; i++) {
         printf ("\n ");
         mpz_out_str (stdout, 10, j [i]);
      }
      printf ("\n");
      printf ("--- Time for j: %.1f\n", cm_timer_get (clock));
   }

   cm_class_clear (&c);

   return j;
}

/*****************************************************************************/
/*****************************************************************************/