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

static void quadratic_basis (mpz_ptr omega, int_cl_t d, mpz_srcptr P);
static void quadraticx_mod_p (mpzx_ptr f, mpzx_srcptr g, mpzx_srcptr h,
   mpz_srcptr omega, mpz_srcptr p);
static void mpzx_eval_mod_p (mpz_ptr val, mpzx_srcptr f, mpz_srcptr x,
   mpz_srcptr p);
static void mpzxx_eval_mod_p (mpzx_ptr val, mpzx_t *f, int deg,
   mpz_srcptr y, mpz_srcptr p);
static void quadraticx_eval_mod_p (mpz_ptr val, mpzx_srcptr g,
   mpzx_srcptr h, mpz_srcptr x, mpz_srcptr omega, mpz_srcptr p);
static void quadraticxx_eval_mod_p (mpzx_ptr val, mpzx_t *g, mpzx_t *h,
   int deg, mpz_srcptr y, mpz_srcptr omega, mpz_srcptr p);
static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose);
static void get_tower_root_mod_p (mpz_ptr root, mpzx_tower_srcptr t,
   mpz_srcptr p, bool verbose);
static void get_quadratic_tower_root_mod_p (mpz_ptr root,
   mpzx_tower_srcptr t, mpzx_tower_srcptr u, mpz_srcptr omega,
   mpz_srcptr p, bool verbose);
static mpz_t* get_j_mod_P_from_modular (int *no, const char* modpoldir,
   char type, int level, mpz_t root, mpz_t P);
static mpz_t* simpleeta_cm_get_j_mod_P (cm_class_t c, mpz_t root, mpz_t P,
   int *no);

/*****************************************************************************/

static void quadratic_basis (mpz_ptr omega, int_cl_t d, mpz_srcptr P)
   /* Compute the second element omega of the standard quadratic basis
      for the fundamental discriminant d modulo P, which is assumed to
      split in the quadratic field. */
{
   cm_nt_mpz_tonelli (omega, d, P);
   if (d % 4 != 0)
      mpz_add_ui (omega, omega, 1ul);
      /* Since sqrt(d)!=p-1 mod p, omega is still <p. */
   if (!mpz_divisible_2exp_p (omega, 1))
      mpz_add (omega, omega, P);
   mpz_divexact_ui (omega, omega, 2);
}

/*****************************************************************************/

static void quadraticx_mod_p (mpzx_ptr f, mpzx_srcptr g, mpzx_srcptr h,
   mpz_srcptr omega, mpz_srcptr p)
   /* Compute f = g + omega * h modulo p, where f is assumed to be
      different from g and h and g and h have the same degree. */
{
   mpz_t tmp;
   int i;

   mpz_init (tmp);

   for (i = 0; i <= g->deg; i++) {
      mpz_mod (f->coeff [i], g->coeff [i], p);
      mpz_mod (tmp, h->coeff [i], p);
      mpz_mul (tmp, tmp, omega);
      mpz_mod (tmp, tmp, p);
      mpz_add (f->coeff [i], f->coeff [i], tmp);
      mpz_mod (f->coeff [i], f->coeff [i], p);
   }

   mpz_clear (tmp);
}

/*****************************************************************************/

static void mpzx_eval_mod_p (mpz_ptr val, mpzx_srcptr f, mpz_srcptr x,
   mpz_srcptr p)
   /* Compute f(x) mod p using a Horner scheme; for better efficiency,
      arg should be an element of [0, p-1]. */
{
   mpz_t res;
   int i;

   mpz_init (res);

   if (f->deg == -1)
      mpz_set_ui (res, 0);
   else {
      mpz_mod (res, f->coeff [f->deg], p);
      for (i = f->deg - 1; i >= 0; i--) {
         mpz_mul (res, res, x);
         mpz_add (res, res, f->coeff [i]);
         mpz_mod (res, res, p);
      }
   }

   mpz_set (val, res);

   mpz_clear (res);
}

/*****************************************************************************/

static void mpzxx_eval_mod_p (mpzx_ptr val, mpzx_t *f, int deg,
   mpz_srcptr y, mpz_srcptr p)
   /* Consider f as a polynomial of degree deg in a variable x with
      coefficients that are polynomials in y. Evaluate f in the variable y
      modulo p. Assume that val is already initialised of degree deg. */
{
   int i;

   for (i = 0; i <= deg; i++)
      mpzx_eval_mod_p (val->coeff [i], f [i], y, p);
}

/*****************************************************************************/

static void quadraticx_eval_mod_p (mpz_ptr val, mpzx_srcptr g,
   mpzx_srcptr h, mpz_srcptr x, mpz_srcptr omega, mpz_srcptr p)
   /* Compute f(x) mod p, where f=g+omega*h. */
{
   mpz_t tmp;

   mpz_init (tmp);

   mpzx_eval_mod_p (val, g, x, p);
   mpzx_eval_mod_p (tmp, h, x, p);
   mpz_mul (tmp, tmp, omega);
   mpz_add (val, val, tmp);
   mpz_mod (val, val, p);

   mpz_clear (tmp);
}

/*****************************************************************************/

static void quadraticxx_eval_mod_p (mpzx_ptr val, mpzx_t *g, mpzx_t *h,
   int deg, mpz_srcptr y, mpz_srcptr omega, mpz_srcptr p)
   /* Consider f=g+omega*h, where g and h are polynomials of degree deg in
      a variable x with coefficients that are polynomials in y. Evaluate f
      in the variable y modulo p. Assume that val is already initialised
      of degree deg. */
{
   int i;

   for (i = 0; i <= deg; i++)
      quadraticx_eval_mod_p (val->coeff [i], g [i], h [i], y, omega, p);
}

/*****************************************************************************/

static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose)
   /* Return a root of the minimal polynomial modulo P in root. */

{
   cm_timer clock;
   mpz_t omega;
   mpzx_t minpoly_P;

   if (c.field == CM_FIELD_REAL)
      cm_pari_oneroot (root, c.minpoly, P, verbose);
   else {
      mpz_init (omega);
      cm_timer_start (clock);
      quadratic_basis (omega, c.dfund, P);
      cm_timer_stop (clock);
      if (verbose)
         printf ("--- Time for square root: %.1f\n", cm_timer_get (clock));

      mpzx_init (minpoly_P, c.minpoly->deg);
      quadraticx_mod_p (minpoly_P, c.minpoly, c.minpoly_complex, omega, P);
      cm_pari_oneroot (root, minpoly_P, P, verbose);

      mpz_clear (omega);
      mpzx_clear (minpoly_P);
   }

   if (verbose) {
      printf ("Root: ");
      mpz_out_str (stdout, 10, root);
      printf ("\n");
   }
}

/*****************************************************************************/

static void get_tower_root_mod_p (mpz_ptr root, mpzx_tower_srcptr t,
   mpz_srcptr p, bool verbose)
   /* Let t be a class field tower in which p splits totally. Return a root
      mod p of the class polynomial behind the tower. */
{
   int i;
   mpzx_t fp;

   cm_pari_oneroot (root, t->W [0][0], p, verbose);
   for (i = 1; i < t->levels; i++) {
      mpzx_init (fp, t->d [i]);
      mpzxx_eval_mod_p (fp, t->W [i], t->d [i], root, p);
      cm_pari_oneroot (root, fp, p, verbose);
      mpzx_clear (fp);
   }

   if (verbose) {
      printf ("Root: ");
      mpz_out_str (stdout, 10, root);
      printf ("\n");
   }
}

/*****************************************************************************/

static void get_quadratic_tower_root_mod_p (mpz_ptr root,
   mpzx_tower_srcptr t, mpzx_tower_srcptr u, mpz_srcptr omega, mpz_srcptr p,
   bool verbose)
   /* Assume that t+omega*u is a class field tower defined over the
      imaginary-quadratic field defined by omega in which p splits
      completely. Return a root modulo p of the class polynomial behind
      the tower. */
{
   int i;
   mpzx_t fp;

   mpzx_init (fp, t->d [0]);
   quadraticx_mod_p (fp, t->W [0][0], u->W [0][0], omega, p);
   cm_pari_oneroot (root, fp, p, verbose);
   mpzx_clear (fp);
   for (i = 1; i < t->levels; i++) {
      mpzx_init (fp, t->d [i]);
      quadraticxx_eval_mod_p (fp, t->W [i], u->W [i], t->d [i], root,
         omega, p);
      cm_pari_oneroot (root, fp, p, verbose);
      mpzx_clear (fp);
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
   mpzx_t poly_j;
      /* a polynomial one of whose roots is j mod P */

   cm_modpol_read_specialised_mod (poly_j, level, type, P, root, modpoldir);
   j = cm_pari_find_roots (no, poly_j, P);
   mpzx_clear (poly_j);

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
   const char* modpoldir, bool pari, bool readwrite, bool tower,
   bool verbose)
   /* pari indicates whether the pari stack should be initialised; if set
      to false, this is the responsability of the calling function.
      If readwrite is true, tries to read the polynomial from a file, and    */
   /* computes it and writes it to the file if it is not yet present.        */

{
   cm_class_t c;
   mpz_t *j;
   mpz_t root, d_mpz, tmp, tmp2, f24;
   cm_timer clock;

   cm_class_init (&c, d, inv, pari, verbose);

   /* In the tower case, reading and writing to files is not yet
      implemented; compute the class polynomial then. */
   if (tower || !readwrite || !cm_class_read (c)) {
      if (tower)
         cm_class_compute_minpoly (c, false, true, readwrite, false, verbose);
      else
         cm_class_compute_minpoly (c, true, false, readwrite, false, verbose);
   }

   cm_timer_start (clock);
   mpz_init (root);
   if (!tower)
      get_root_mod_P (c, root, P, verbose);
   else {
      if (c.field == CM_FIELD_REAL)
         get_tower_root_mod_p (root, c.tower, P, verbose);
      else {
         mpz_t omega;
         mpz_init (omega);
         quadratic_basis (omega, c.dfund, P);
         get_quadratic_tower_root_mod_p (root, c.tower, c.tower_complex,
            omega, P, verbose);
         mpz_clear (omega);
      }
   }

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

         if (d % 3 == 0)
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
