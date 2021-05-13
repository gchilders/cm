/*

class.c - code for computing class polynomials

Copyright (C) 2009, 2010, 2011, 2012, 2015, 2016, 2017, 2018, 2021 Andreas Enge

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

static double class_get_valuation (cm_class_srcptr c);
   /* return some value related to heights and depending on the function     */
static int class_get_height (cm_class_srcptr c);
   /* in the real case, returns the binary length of the largest             */
   /* coefficient of the minimal polynomial                                  */
   /* in the complex case, returns the binary length of the largest          */
   /* coefficient with respect to the decomposition over an integral basis   */

static void correct_nsystem_entry (cm_form_t *Q, int_cl_t N, int_cl_t b0,
   int_cl_t d);
   /* changes Q to be compatible with the N-system condition */
static void compute_nsystem (cm_form_t *nsystem, int *conj, cm_class_srcptr c,
   bool verbose);
   /* computes and returns an N-system of forms, together with information
      on the embeddings in the case of real class polynomials in conj */
static fprec_t compute_precision (cm_class_srcptr c, bool verbose);

static void eval (cm_class_srcptr c, cm_modclass_t mc, ctype rop, cm_form_t Q);
static void compute_conjugates (ctype *conjugate, cm_form_t *nsystem,
   int *conj, cm_class_srcptr c, cm_modclass_t mc, bool verbose);


/*****************************************************************************/
/*                                                                           */
/* constructor and destructor                                                */
/*                                                                           */
/*****************************************************************************/

void cm_class_init (cm_class_ptr c, cm_param_srcptr param, bool pari,
   bool verbose)

{
   int i;
   int one [] = {1};

   c->invariant = param->invariant;
   c->dfund = cm_classgroup_fundamental_discriminant (param->d);
   for (i = 0; i < 6; i++)
      c->p [i] = param->p [i];
   c->e = param->e;
   c->s = param->s;
   strncpy (c->paramstr, param->str, 255);
   if (verbose)
      printf ("\nDiscriminant %"PRIicl", invariant %c, parameter %s\n",
               param->d, c->invariant, c->paramstr);

   c->field = param->field;

   c->pari = pari;
   if (pari) {
      pari_init_opts (1000000, 0, INIT_JMPm | INIT_DFTm);
         /* Do not capture SIGSEGV. */
      paristack_setsize (1000000, 1000000000);
   }

   cm_classgroup_init (&(c->cl), param->d, verbose);
   mpzx_init (c->minpoly, c->cl.h);
   if (c->field == CM_FIELD_COMPLEX)
      mpzx_init (c->minpoly_complex, c->cl.h);

   if (c->cl.h == 1) {
      mpzx_tower_init (c->tower, 1, one);
      if (c->field == CM_FIELD_COMPLEX)
         mpzx_tower_init (c->tower_complex, 1, one);
   }
   else {
      mpzx_tower_init (c->tower, c->cl.levels, c->cl.deg);
      if (c->field == CM_FIELD_COMPLEX)
         mpzx_tower_init (c->tower_complex, c->cl.levels, c->cl.deg);
   }
}

/*****************************************************************************/

void cm_class_clear (cm_class_ptr c)

{
   mpzx_clear (c->minpoly);
   mpzx_tower_clear (c->tower);
   if (c->field == CM_FIELD_COMPLEX) {
      mpzx_tower_clear (c->tower_complex);
   }

   if (c->pari)
      pari_close ();

   cm_classgroup_clear (&(c->cl));

   ffree_cache ();
}


/*****************************************************************************/
/*                                                                           */
/* valuation at infinity and height                                          */
/*                                                                           */
/*****************************************************************************/

static double class_get_valuation (cm_class_srcptr c)
   /* returns the (negative) order of the modular function at infinity       */

{
   double result;

   switch (c->invariant) {
   case CM_INVARIANT_J:
      result = 1;
      break;
   case CM_INVARIANT_GAMMA2:
      result = 1.0 / 3;
      break;
   case CM_INVARIANT_GAMMA3:
      result = 0.5;
      break;
   case CM_INVARIANT_ATKIN:
      if (c->p [0] == 47)
         result = 1.0 / 24;
      else if (c->p [0] == 59)
         result = 1.0 / 30;
      else if (c->p [0] == 71)
         result = 1.0 / 36;
      else /* 131 */
         result = 1.0 / 33;
      break;
   case CM_INVARIANT_WEBER:
      result = 1.0 / 72;
      break;
   case CM_INVARIANT_DOUBLEETA:
   case CM_INVARIANT_MULTIETA:
   {
      int num = 1, den = 1, i;
      for (i = 0; c->p [i] != 0; i++) {
         num *= c->p [i] - 1;
         den *= c->p [i] + 1;
      }
      if (i == 2)
         result = num / (double) (12 * den);
      else if (i == 3)
         result = num / (double) (6 * den);
      else /* i == 4 */
         result = num / (double) (3 * den);
   }
      break;
   case CM_INVARIANT_SIMPLEETA:
      result = (c->p [0] - 1) / (double) (24 * (c->p [0] + 1));
      break;
   default: /* should not occur */
      printf ("class_get_valuation called for unknown class ");
      printf ("invariant\n");
      exit (1);
   }

   result *= c->e;

   return result;
}

/*****************************************************************************/

static int class_get_height (cm_class_srcptr c)
   /* In the real case, return the binary length of the largest coefficient
      of the minimal polynomial; in the complex case, return the binary
      length of the largest coefficient with respect to the decomposition
      over an integral basis of the imaginary-quadratic field. */
{
   int i, height, cand;

   height = -1;
   for (i = 0; i < c->minpoly->deg; i++) {
      cand = mpz_sizeinbase (c->minpoly->coeff [i], 2);
      if (cand > height)
         height = cand;
   }
   if (c->field == CM_FIELD_COMPLEX)
      for (i = 0; i < c->minpoly_complex->deg; i++) {
         cand = mpz_sizeinbase (c->minpoly_complex->coeff [i], 2);
         if (cand > height)
            height = cand;
   }

   return height;
}

/*****************************************************************************/
/*                                                                           */
/* computing the class polynomial                                            */
/*                                                                           */
/*****************************************************************************/

static void correct_nsystem_entry (cm_form_t *Q, int_cl_t N, int_cl_t b0,
   int_cl_t d)
   /* Changes the form Q by a unimodular transformation so that Q.a is
      coprime to N and Q.b is congruent to b0 modulo 2*N. */

{
   int_cl_t c, tmp;

   /* First achieve gcd (Q->a, N) = 1, which is likely to hold already.   */
   c = (Q->b * Q->b - d) / (4 * Q->a) ;
   if (cm_classgroup_gcd (Q->a, N) != 1) {
      /* Translation by k yields C' = A k^2 + B k + C; we wish to reach   */
      /* gcd (C', N) = 1, so for each prime p dividing N, this excludes   */
      /* at most two values modulo p. For p = 2, A and B odd and C even,  */
      /* there is no solution; in this case, we first apply S to exchange */
      /* A and C.                                                         */
      if (N % 2 == 0 && Q->a % 2 != 0 && Q->b % 2 != 0 && c % 2 == 0) {
         tmp = Q->a;
         Q->a = c;
         c = tmp;
         Q->b = -Q->b;
      }
      while (cm_classgroup_gcd (c, N) != 1) {
         /* Translate by 1 */
         c += Q->a + Q->b;
         Q->b += 2 * Q->a;
      }
      /* Apply S */
      tmp = Q->a;
      Q->a = c;
      c = tmp;
      Q->b = -Q->b;
   }
   /* Translate so that Q->b = b0 mod (2 N).                              */
   while ((Q->b - b0) % (2*N) != 0) {
      c += Q->a + Q->b;
      Q->b += 2 * Q->a;
   }
}

/*****************************************************************************/

static void compute_nsystem (cm_form_t *nsystem, int *conj, cm_class_srcptr c,
   bool verbose)
   /* Compute an N-system, or to be more precise, some part of an N-system
      that yields all different conjugates up to complex conjugation;
      this information is passed in conj as follows:
      If the class polynomial is real, then conj [i] == j if form j in
      If the class polynomial is complex, then conj [i] = i. */

{
   int_cl_t b0, N;
   cm_form_t neutral, inverse;
   int h1, h2;
   int i, j;

   /* Compute the targeted b0 for the N-system and (in the real case) the
      neutral form. */
   if (c->invariant == CM_INVARIANT_SIMPLEETA) {
      bool ok = false;

      if (c->p[0] != 4) {
         b0 = c->cl.d % 2;
         while (!ok) {
            b0 += 2;
            if ((b0*b0 - c->cl.d) % (4*c->p[0]) != 0)
               ok = false;
            else if (c->p[0] != 2 && (c->s/c->e) % 2 == 0 && (b0 - 1) % 4 != 0)
               ok = false;
            else if (c->p[0] != 2 && (c->s/c->e) % 3 == 0 && b0 % 3 != 0)
               ok = false;
            else
               ok = true;
         }
      }
      else {
        b0 = (int_cl_t) -7;
      }
      if (verbose)
         printf ("N %i\ns %i\ne %i\nb0 %"PRIicl"\n", c->p[0], c->s, c->e, b0);
      N = c->p[0] * c->s / c->e;
   }

   else {
      neutral.a = 1;
      if (c->cl.d % 2 == 0)
         neutral.b = 0;
      else
         neutral.b = 1;

      switch (c->invariant) {
         case CM_INVARIANT_J:
            b0 = c->cl.d % 2;
            /* An even N makes c even if 2 is split so that during the
               evaluation of eta(z/2) for f1 all forms can be taken from
               the precomputed ones. */
            N = 2;
            break;
         case CM_INVARIANT_GAMMA2:
            b0 = 3 * (c->cl.d % 2);
            /* Use an even N as for j. */
            N = 6;
            break;
         case CM_INVARIANT_GAMMA3:
            b0 = 1;
            N = 2;
            break;
         case CM_INVARIANT_ATKIN:
            N = c->p [0];
            if (c->cl.d % 2 == 0)
               b0 = 0;
            else
               b0 = 1;
            while ((b0*b0 - c->cl.d) % N != 0)
               b0 += 2;
            neutral.a = N;
            neutral.b = -b0;
            cm_classgroup_reduce (&neutral, c->cl.d);
            break;
         case CM_INVARIANT_WEBER:
            neutral.a = 1;
            neutral.b = 0;
            b0 = 0;
            N = 48;
            break;
         case CM_INVARIANT_DOUBLEETA:
         case CM_INVARIANT_MULTIETA:
         {
            int_cl_t C;
            N = 1;
            for (i = 0; c->p [i] != 0; i++)
               N *= c->p [i];
            if (c->cl.d % 2 == 0)
               b0 = 2;
            else
               b0 = 1;
            while (true) {
               C = (b0*b0 - c->cl.d) / 4;
               if (C % N == 0 && cm_nt_gcd (C / N, N) == 1)
                  break;
               b0 += 2;
            }
            neutral.a = N;
            neutral.b = -b0;
            cm_classgroup_reduce (&neutral, c->cl.d);
            /* The neutral form corresponds to the product of the primes,    */
            /* but the n-system needs to take s/e into account               */
            N *= c->s / c->e;
            break;
         }
         default: /* should not occur */
            printf ("compute_nsystem called for unknown class invariant\n");
            exit (1);
      }
   }

   for (i = 0; i < c->cl.h; i++) {
      nsystem [i] = c->cl.form [i];
      conj [i] = -1;
   }

   for (i = 0; i < c->cl.h; i++)
      /* Pair forms yielding complex conjugate roots. */
      if (c->field == CM_FIELD_REAL) {
         if (conj [i] == -1) {
            /* form did not yet occur in a pair; look for its inverse
               with respect to neutral_class */
            nsystem [i].b = -nsystem [i].b;
            cm_classgroup_compose (&inverse, neutral, nsystem [i], c->cl.d);
            nsystem [i].b = -nsystem [i].b;
            j = 0;
            /* So far, nsystem still contains the reduced forms, so we may look
               for the inverse form. */
            while (nsystem [j].a != inverse.a || nsystem [j].b != inverse.b)
               j++;
            conj [i] = j;
            conj [j] = i;
         }
      }
      else /* c->field == CM_FIELD_COMPLEX */
         conj [i] = i;

   /* Now modify the entries of nsystem. */
   for (i = 0; i < c->cl.h; i++)
      if (conj [i] >= i)
         correct_nsystem_entry (&(nsystem [i]), N, b0, c->cl.d);

   /* Compute h1 and h2, only for printing. */
   if (verbose) {
      h1 = 0;
      h2 = 0;
      for (i = 0; i < c->cl.h; i++)
         if (conj [i] == i)
            h1++;
         else if (conj [i] > i)
            h2++;
      printf ("h = %i, h1 = %i, h2 = %i\n", c->cl.h, h1, h2);
   }
}

/*****************************************************************************/

static fprec_t compute_precision (cm_class_srcptr c, bool verbose) {
   /* Returns an approximation of the required precision. */

   const double C = 2114.567;
   const double pisqrtd
      = 3.14159265358979323846 * sqrt ((double) (-c->cl.d));
   const double cf = class_get_valuation (c);
   double x, binom = 1.0, prec = 0, M;
   int_cl_t amax;
   int i, m;
   fprec_t precision;

   /* formula of Lemma 8 of [Sutherland11] */
   amax = 0;
   for (i = 0; i < c->cl.h; i++) {
      x = pisqrtd / c->cl.form [i].a;
      if (x < 42)
         M = log (exp (x) + C);
      else /* prevent overflow in exponential without changing the result */
         M = x;
      prec += M;
      if (c->cl.form [i].a > amax)
         amax = c->cl.form [i].a;
   }
   M = exp (pisqrtd / amax) + C;
   m = (int) ((c->cl.h + 1) / (M + 1));
   for (i = 1; i <= m; i++)
      binom *= (double) (c->cl.h - 1 + i) / i / M;
   prec = ceil ((prec + log (binom)) / log (2.0) * cf);

   if (c->invariant == CM_INVARIANT_GAMMA3) {
      /* Increase the height estimate by the bit size of sqrt (|D|)^h in
         the constant coefficient.*/
      prec += (int) (log ((double) (-c->cl.d)) / log (2.0) / 2.0 * c->cl.h);
      if (verbose)
         printf ("Corrected bound for gamma3:     %ld\n", (long int) prec);
   }

   /* add a security margin */
   precision = (fprec_t) (prec + 256);

   if (verbose)
      printf ("Precision:                      %ld\n", (long int) precision);

   return precision;
}
/*****************************************************************************/

static void eval (cm_class_srcptr c, cm_modclass_t mc, ctype rop, cm_form_t Q)

{
   switch (c->invariant) {
   case CM_INVARIANT_J:
      cm_modclass_j_eval_quad (mc, rop, Q.a, Q.b);
      break;
   case CM_INVARIANT_GAMMA2:
      cm_modclass_gamma2_eval_quad (mc, rop, Q.a, Q.b);
      break;
   case CM_INVARIANT_GAMMA3:
      cm_modclass_gamma3_eval_quad (mc, rop, Q.a, Q.b);
      break;
   case CM_INVARIANT_ATKIN:
      cm_modclass_atkinhecke_level_eval_quad (mc, rop, Q.a, Q.b, c->p [0]);
      break;
   case CM_INVARIANT_SIMPLEETA:
   case CM_INVARIANT_DOUBLEETA:
   case CM_INVARIANT_MULTIETA:
         cm_modclass_multieta_eval_quad (mc, rop, Q.a, Q.b, c->p, c->e);
      break;
   case CM_INVARIANT_WEBER:
      if (c->p [0] == 1) {
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 2);
         cmul_fr (rop, rop, mc.sqrt2_over2);
      }
      else if (c->p [0] == 3)
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 1);
      else if (c->p [0] == 5) {
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 2);
         csqr (rop, rop);
         cdiv_ui (rop, rop, 2ul);
      }
      else if (c->p [0] == 7) {
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 1);
         cmul_fr (rop, rop, mc.sqrt2_over2);
      }
      else if (c->p [0] == 2 || c->p [0] == 6) {
         cm_modclass_f1_eval_quad (mc, rop, Q.a, Q.b, 1);
         csqr (rop, rop);
         cmul_fr (rop, rop, mc.sqrt2_over2);
      }
      else {
         /* c->p [0] == 4 */
         cm_modclass_f1_eval_quad (mc, rop, Q.a, Q.b, 1);
         cpow_ui (rop, rop, 4ul);
         cmul_fr (rop, rop, mc.sqrt2_over4);
      }

      if (c->cl.d % 3 == 0)
         cpow_ui (rop, rop, 3ul);

      if (c->p [0] != 3 && c->p [0] != 5)
         if (cm_classgroup_kronecker ((int_cl_t) 2, Q.a) == -1)
            cneg (rop, rop);

      break;
   default: /* should not occur */
      printf ("class_eval called for unknown class invariant\n");
      exit (1);
   }
}

/*****************************************************************************/

static void compute_conjugates (ctype *conjugate, cm_form_t *nsystem,
   int *conj, cm_class_srcptr c, cm_modclass_t mc, bool verbose)
   /* computes the conjugates of the singular value over Q */

{
   int i;

   for (i = 0; i < c->cl.h; i++) {
      if (conj [i] >= i)
         eval (c, mc, conjugate [i], nsystem [i]);
      if (verbose && i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   if (verbose)
      printf ("\n");
}

/*****************************************************************************/

bool cm_class_compute_minpoly (cm_class_ptr c, bool classpol, bool tower,
   bool disk, bool print, bool verbose)
   /* At least one of classpol and tower needs to be set to true:
      classpol indicates whether the (absolute) class polynomial should be
      computed; tower indicates whether the class polynomial should be
      decomposed as a Galois tower.
      disk indicates whether the result should be written to disk.
      print indicates whether the result should be printed on screen.
      The return value reflects the success of the computation. */
{
   cm_form_t *nsystem;
   int *conj;
   fprec_t prec;
   cm_modclass_t mc;
   ctype *conjugate;
   mpfrx_t mpol;
   mpcx_t mpolc;
   mpfrx_tower_t t;
   mpcx_tower_t tc;
   int i;
   bool ok = true;
   cm_timer clock_global, clock_local;

   if (!classpol && !tower) {
      printf ("***** Error: cm_class_compute_minpoly called with nothing "
              "to compute\n");
      return false;
   }

   cm_timer_start (clock_global);

   nsystem = (cm_form_t *) malloc (c->cl.h * sizeof (cm_form_t));
   conj = (int *) malloc (c->cl.h * sizeof (int));
   cm_timer_start (clock_local);
   compute_nsystem (nsystem, conj, c, verbose);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for N-system: %.1f\n", cm_timer_get (clock_local));
   prec = compute_precision (c, verbose);

   conjugate = (ctype *) malloc (c->cl.h * sizeof (ctype));
   for (i = 0; i < c->cl.h; i++)
      if (conj [i] >= i)
         cinit (conjugate [i], prec);
   cm_timer_start (clock_local);
   cm_modclass_init (&mc, c->cl, prec, verbose);
   compute_conjugates (conjugate, nsystem, conj, c, mc, verbose);
   cm_modclass_clear (&mc);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for conjugates: %.1f\n", cm_timer_get (clock_local));

   if (classpol) {
      cm_timer_start (clock_local);
      if (c->field == CM_FIELD_REAL) {
         mpfrx_init (mpol, c->minpoly->deg + 1, prec);
         mpfrcx_reconstruct_from_roots (mpol, conjugate, conj,
            c->minpoly->deg);
         ok &= cm_mpfrx_get_mpzx (c->minpoly, mpol);
         if (print && ok) {
            mpzx_print_pari (stdout, c->minpoly, NULL);
            printf ("\n");
         }
         mpfrx_clear (mpol);
      }
      else {
         mpcx_init (mpolc, c->minpoly->deg + 1, prec);
         mpcx_reconstruct_from_roots (mpolc, conjugate, c->minpoly->deg);
         ok &= cm_mpcx_get_quadraticx (c->minpoly, c->minpoly_complex,
            mpolc, c->dfund);
         mpcx_clear (mpolc);
         if (print && ok) {
            printf ("(");
            mpzx_print_pari (stdout, c->minpoly, NULL);
            printf (")+o*(");
            mpzx_print_pari (stdout, c->minpoly_complex, NULL);
            printf (")\n");
         }
      }
      cm_timer_stop (clock_local);
      if (verbose)
         printf ("--- Time for minimal polynomial reconstruction: %.1f\n",
                 cm_timer_get (clock_local));
   }

   if (tower && ok) {
      cm_timer_start (clock_local);
      if (c->field == CM_FIELD_REAL) {
         mpfrx_tower_init (t, c->tower->levels, c->tower->d, prec);
         mpfrcx_tower_decomposition (t, conjugate, conj);
         /* There should be a possibility to save the field tower on disk. */
         ok &= cm_mpfrx_tower_get_mpzx_tower (c->tower, t);
         if (print && ok)
            mpzx_tower_print_pari (stdout, c->tower, "f", NULL);
         mpfrx_tower_clear (t);
      }
      else {
         mpcx_tower_init (tc, c->tower->levels, c->tower->d, prec);
         mpcx_tower_decomposition (tc, conjugate);
         ok = cm_mpcx_tower_get_quadratic_tower (c->tower, c->tower_complex,
            tc, c->dfund);
         mpcx_tower_clear (tc);
         if (print && ok) {
            mpzx_tower_print_pari (stdout, c->tower, "f", NULL);
            mpzx_tower_print_pari (stdout, c->tower_complex, "g", NULL);
         }
      }
      cm_timer_stop (clock_local);
      if (verbose)
         printf ("--- Time for field tower decomposition: %.1f\n",
                 cm_timer_get (clock_local));
   }

   for (i = 0; i < c->cl.h; i++)
      if (conj [i] >= i)
         cclear (conjugate [i]);
   free (conjugate);
   free (nsystem);
   free (conj);

   cm_timer_stop (clock_global);
   if (verbose) {
      printf ("--- Total time for minimal polynomial: %.1f\n",
            cm_timer_get (clock_global));
      if (classpol)
         printf ("Height of minimal polynomial: %d\n",
            class_get_height (c));
   }
   if (disk && ok)
      ok &= cm_class_write (c);

   return ok;
}

/*****************************************************************************/
/*****************************************************************************/
