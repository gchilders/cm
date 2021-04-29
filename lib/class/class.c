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

#include "cm_class-impl.h"

static double class_get_valuation (cm_class_t c);
   /* return some value related to heights and depending on the function     */
static int class_get_height (cm_class_t c);
   /* in the real case, returns the binary length of the largest             */
   /* coefficient of the minimal polynomial                                  */
   /* in the complex case, returns the binary length of the largest          */
   /* coefficient with respect to the decomposition over an integral basis   */

static bool cm_class_compute_parameter (cm_class_t *c, bool verbose);
static void correct_nsystem_entry (cm_form_t *Q, int_cl_t N, int_cl_t b0,
   int_cl_t d);
   /* changes Q to be compatible with the N-system condition */
static void compute_nsystem (cm_form_t *nsystem, int *conj, cm_class_t *c,
   cm_classgroup_t cl, bool verbose);
   /* computes and returns an N-system of forms, together with information
      on the embeddings in the case of real class polynomials in conj */
static fprec_t compute_precision (cm_class_t c, cm_classgroup_t cl,
   bool verbose);

static void eval (cm_class_t c, cm_modclass_t mc, ctype rop, cm_form_t Q);
static void compute_conjugates (ctype *conjugate, cm_form_t *nsystem,
   int *conj, cm_class_t c, cm_modclass_t mc, bool verbose);
static void write_conjugates (cm_class_t c, ctype *conjugates, int *conj);
static bool read_conjugates (cm_class_t c, ctype *conjugates, int *conj);

static bool get_quadratic (mpz_t out1, mpz_t out2, ctype in, int_cl_t d);

static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose);
static mpz_t* cm_get_j_mod_P_from_modular (int *no, const char* modpoldir,
   char type, int level, mpz_t root, mpz_t P);

/* the remaining functions are put separately as their code is quite long    */
static void real_compute_minpoly (cm_class_t c, ctype *conjugate, int *conj,
   bool print);
static void complex_compute_minpoly (cm_class_t c, ctype *conjugate,
   bool print);
static bool doubleeta_compute_parameter (cm_class_t *c);
static mpz_t* weber_cm_get_j_mod_P (cm_class_t c, mpz_t root, mpz_t P,
   int *no, bool verbose);
static mpz_t* simpleeta_cm_get_j_mod_P (cm_class_t c, mpz_t root, mpz_t P,
   int *no);

/*****************************************************************************/
/*                                                                           */
/* constructor and destructor                                                */
/*                                                                           */
/*****************************************************************************/

void cm_class_init (cm_class_t *c, int_cl_t d, char inv, bool verbose)

{
   cm_classgroup_t cl;
   int deg;
   int i, j, k;

   if (d >= 0) {
      printf ("\n*** The discriminant must be negative.\n");
      exit (1);
   }
   else if (d % 4 != 0 && (d - 1) % 4 != 0) {
      printf ("\n*** %"PRIicl" is not a quadratic discriminant.\n", d);
      exit (1);
   }
   else {
      c->invariant = inv;
      c->d = d;
      if (!cm_class_compute_parameter (c, verbose))
         exit (1);
      if (inv == CM_INVARIANT_WEBER && c->p [1] == 1)
         /* special case Weber with d=1 (4): compute ring class field for 4d */
         /* If d=1 (8), this is the same as the Hilbert class field;         */
         /* if d=5 (8), it is a degree 3 relative extension, which will be   */
         /* corrected when computing a CM curve by applying a 2-isogeny/     */
         c->d = 4 * d;
   }
   if (verbose)
      printf ("\nDiscriminant %"PRIicl", invariant %c, parameter %s\n",
               c->d, c->invariant, c->paramstr);

   if (inv == CM_INVARIANT_SIMPLEETA)
      c->field = CM_FIELD_COMPLEX;
   else if (inv == CM_INVARIANT_MULTIETA) {
      for (i = 0; c->p [i] != 0; i++);
      if (i % 2 == 0)
         c->field = CM_FIELD_REAL;
      else
         c->field = CM_FIELD_COMPLEX;
   }
   else
      c->field = CM_FIELD_REAL;

   c->h = cm_classgroup_h (c->d);
   c->minpoly_deg = c->h;
   c->minpoly = (mpz_t *) malloc ((c->minpoly_deg + 1) * sizeof (mpz_t));
   for (i = 0; i <= c->minpoly_deg; i++)
      mpz_init (c->minpoly [i]);
   if (c->field == CM_FIELD_COMPLEX) {
      c->minpoly_complex = (mpz_t *) malloc (c->minpoly_deg * sizeof (mpz_t));
      for (i = 0; i < c->minpoly_deg; i++)
         mpz_init (c->minpoly_complex [i]);
   }

   cm_classgroup_init (&cl, c->d, false);
   if (cl.levels == 0) {
      deg = 1;
      c->tower_levels = 1;
   }
   else {
      deg = cl.deg [0];
      c->tower_levels = cl.levels;
   }
   c->tower_d = (int *) malloc (c->tower_levels * sizeof (int));
   c->tower = (mpz_t ***) malloc (c->tower_levels * sizeof (mpz_t **));
   c->tower_d [0] = deg;
   c->tower [0] = (mpz_t **) malloc (1 * sizeof (mpz_t *));
   c->tower [0][0] = (mpz_t *) malloc ((deg + 1) * sizeof (mpz_t));
   for (k = 0; k <= deg; k++)
      mpz_init (c->tower [0][0][k]);
   for (i = 1; i < c->tower_levels; i++) {
      c->tower_d [i] = cl.deg [i];
      c->tower [i] = (mpz_t **) malloc ((cl.deg [i] + 1) * sizeof (mpz_t *));
      for (j = 0; j <= cl.deg [i]; j++) {
         c->tower [i][j] = (mpz_t *) malloc (deg * sizeof (mpz_t));
         for (k = 0; k < deg; k++)
            mpz_init (c->tower [i][j][k]);
      }
      deg *= cl.deg [i];
   }
   if (c->field == CM_FIELD_COMPLEX) {
      c->tower_complex = (mpz_t ***) malloc (c->tower_levels * sizeof (mpz_t **));
      if (cl.levels == 0)
         deg = 1;
      else
         deg = cl.deg [0];
      c->tower_complex [0] = (mpz_t **) malloc (1 * sizeof (mpz_t *));
      c->tower_complex [0][0]
         = (mpz_t *) malloc ((deg + 1) * sizeof (mpz_t));
      for (k = 0; k <= deg; k++)
         mpz_init (c->tower_complex [0][0][k]);
      for (i = 1; i < cl.levels; i++) {
         c->tower_complex [i]
            = (mpz_t **) malloc ((cl.deg [i] + 1) * sizeof (mpz_t *));
         for (j = 0; j <= cl.deg [i]; j++) {
            c->tower_complex [i][j]
               = (mpz_t *) malloc (deg * sizeof (mpz_t));
            for (k = 0; k < deg; k++)
               mpz_init (c->tower_complex [i][j][k]);
         }
         deg *= cl.deg [i];
      }
   }
   cm_classgroup_clear (&cl);
}

/*****************************************************************************/

void cm_class_clear (cm_class_t *c)

{
   int d;
   int i, j, k;

   for (i = 0; i <= c->minpoly_deg; i++)
      mpz_clear (c->minpoly [i]);
   free (c->minpoly);
   if (c->field == CM_FIELD_COMPLEX) {
      for (i = 0; i < c->minpoly_deg; i++)
         mpz_clear (c->minpoly_complex [i]);
      free (c->minpoly_complex);
   }

   if (c->field == CM_FIELD_COMPLEX) {
      d = c->tower_d [0];
      for (k = 0; k <= d; k++)
         mpz_clear (c->tower_complex [0][0][k]);
      free (c->tower_complex [0][0]);
      free (c->tower_complex [0]);
      for (i = 1; i < c->tower_levels; i++) {
         for (j = 0; j <= c->tower_d [i]; j++) {
            for (k = 0; k < d; k++)
               mpz_clear (c->tower_complex [i][j][k]);
            free (c->tower_complex [i][j]);
         }
         d *= c->tower_d [i];
            /* total degree of the level over \Q */
         free (c->tower_complex [i]);
      }
      free (c->tower_complex);
   }

   d = c->tower_d [0];
   for (k = 0; k <= d; k++)
      mpz_clear (c->tower [0][0][k]);
   free (c->tower [0][0]);
   free (c->tower [0]);
   for (i = 1; i < c->tower_levels; i++) {
      for (j = 0; j <= c->tower_d [i]; j++) {
         for (k = 0; k < d; k++)
            mpz_clear (c->tower [i][j][k]);
         free (c->tower [i][j]);
      }
      d *= c->tower_d [i];
         /* total degree of the level over \Q */
      free (c->tower [i]);
   }
   free (c->tower);
   free (c->tower_d);

   ffree_cache ();
}

/*****************************************************************************/

static bool cm_class_compute_parameter (cm_class_t *c, bool verbose)
      /* tests whether the discriminant is suited for the chosen invariant and  */
      /* in this case computes and returns the parameter p                      */
      /* otherwise, returns 0                                                   */

{
   int i;
   char* pointer;

   switch (c->invariant) {
      case CM_INVARIANT_J:
         c->p [0] = 1;
         c->p [1] = 0;
         c->s = 1;
         c->e = 1;
         break;
      case CM_INVARIANT_GAMMA2:
         if (c->d % 3 == 0) {
            if (verbose) {
               printf ("\n*** %"PRIicl" is divisible by 3, so that gamma2 ",
                        c->d);
               printf ("cannot be used.\n");
            }
            return false;
         }
         c->p [0] = 1;
         c->p [1] = 0;
         c->s = 1;
         c->e = 1;
         break;
      case CM_INVARIANT_GAMMA3:
         if (c->d % 2 == 0) {
            if (verbose) {
               printf ("\n*** %"PRIicl" is divisible by 4, so that gamma3 ",
                        c->d);
               printf ("cannot be used.\n");
            }
            return false;
         }
         c->p [0] = 1;
         c->p [1] = 0;
         c->s = 1;
         c->e = 1;
         break;
      case CM_INVARIANT_WEBER:
         if (c->d % 32 == 0) {
            if (verbose) {
               printf ("\n*** %"PRIicl" is divisible by 32, so that the Weber ",
                        c->d);
               printf ("functions cannot be used.\n");
            }
            return false;
         }

         /* If disc is even, let m = -disc and p [0] = m % 8; if disc is  */
         /* odd, compute p [0] for 4*disc and let p [1] = 1.              */
         if (c->d % 4 == 0) {
            c->p [0] = ((-(c->d)) / 4) % 8;
            c->p [1] = 0;
         }
         else {
            c->p [0] = (-(c->d)) % 8;
            c->p [1] = 1;
            c->p [2] = 0;
         }
         c->s = 24;
         c->e = 1;
         if (c->p [0] == 1 || c->p [0] == 2 || c->p [0] == 6)
            c->e *= 2;
         else if (c->p [0] == 4 || c->p [0] == 5)
            c->e *= 4;
         if (c->d % 3 == 0)
            c->e *= 3;
         break;
      case CM_INVARIANT_ATKIN:
         if (cm_classgroup_kronecker (c->d, (int_cl_t) 71) != -1)
            /* factor 36, T_5 + T_29 + 1 */
            c->p [0] = 71;
         else if (cm_classgroup_kronecker (c->d, (int_cl_t) 131) != -1)
            /* factor 33, T_61 + 1 */
            c->p [0] = 131;
         else if (cm_classgroup_kronecker (c->d, (int_cl_t) 59) != -1)
            /* factor 30, T_5 + T_29 */
            c->p [0] = 59;
         else if (cm_classgroup_kronecker (c->d, (int_cl_t) 47) != -1)
            /* factor 24, -T_17 */
            c->p [0] = 47;
         else {
            if (verbose) {
               printf ("\n*** 47, 59, 71 and 131 are inert for %"PRIicl, c->d);
               printf (", so that atkin cannot be used.\n");
            }
            return false;
         }
         c->p [1] = 0;
         c->s = 1;
         c->e = 1;
         break;
      case CM_INVARIANT_MULTIETA:
         c->p [0] = 3;
         c->p [1] = 5;
         c->p [2] = 7;
         c->p [3] = 0;
         c->s = 1;
         c->e = 1;
         break;
      case CM_INVARIANT_DOUBLEETA:
         if (!doubleeta_compute_parameter (c))
            return false;
         break;
      case CM_INVARIANT_SIMPLEETA:
         c->p [0] = 3;
         if (cm_classgroup_kronecker (c->d, (int_cl_t) (c->p [0])) == -1) {
            if (verbose)
               printf ("*** Unsuited discriminant\n\n");
            return false;
         }

         c->p [1] = 0;
         c->s = 12;
         c->e = 12;
         break;
      default: /* should not occur */
         printf ("class_compute_parameter called for "
                  "unknown class invariant '%c'\n", c->invariant);
         exit (1);
   }

   /* create parameter string */
   pointer = c->paramstr;
   i = 0;
   do {
      pointer += sprintf (pointer, "%i_", c->p [i]);
      i++;
   } while (c->p [i] != 0);
   sprintf (pointer, "%i_%i", c->e, c->s);
   return true;
}

/*****************************************************************************/

static bool doubleeta_compute_parameter (cm_class_t *c)
   /* Compute p1 <= p2 prime following Cor. 3.1 of [EnSc04], that is,        */
   /* - 24 | (p1-1)(p2-1).                                                   */
   /* - p1, p2 are not inert;                                                */
   /* - if p1!=p2, then p1, p2 do not divide the conductor;                  */
   /* - if p1=p2=p, then either p splits or divides the conductor.           */
   /* The case p1=p2=2 of Cor. 3.1 is excluded by divisibility by 24.        */
   /* Minimise with respect to the height factor gained, which is            */
   /* 12 psi (p1*p2) / (p1-1)(p2-1);                                         */
   /* then p1, p2 <= the smallest split prime which is 1 (mod 24).           */
   /* Let p = 1000*p2 + p1.                                                */

{
   int_cl_t cond2 = c->d / cm_classgroup_fundamental_discriminant (c->d);
      /* square of conductor */
   const unsigned long int maxprime = 997;
   unsigned long int primelist [169];
      /* list of suitable primes, terminated by 0; big enough to hold all    */
      /* primes <= maxprime                                                  */
   unsigned long int p1, p2, p1opt = 0, p2opt = 0;
   double quality, opt;
   bool ok;
   int i, j;

   /* determine all non-inert primes */
   i = 0;
   p1 = 2;
   ok = false;
   do {
      int kro = cm_classgroup_kronecker (c->d, (int_cl_t) p1);
      if (kro != -1) {
         primelist [i] = p1;
         i++;
      }
      if (kro == 1 && (p1 - 1) % 24 == 0)
         ok = true;
      else
         p1 = cm_nt_next_prime (p1);
   }
   while (p1 <= maxprime && !ok);
   primelist [i] = 0;

   /* search for the best tuple */
   opt = 0.0;
   for (j = 0, p2 = primelist [j]; p2 != 0; j++, p2 = primelist [j])
      for (i = 0, p1 = primelist [i]; i <= j; i++, p1 = primelist [i])
         if (   ((p1 - 1)*(p2 - 1)) % 24 == 0
             && (   (p1 != p2 && cond2 % p1 != 0 && cond2 % p2 != 0)
                 || (p1 == p2 && ((-c->d) % p1 != 0 || cond2 % p1 == 0)))) {
            quality = (p1 == p2 ? p1 : p1 + 1) * (p2 + 1) / (double) (p1 - 1)
               / (double) (p2 - 1);
            if (quality > opt) {
               p1opt = p1;
               p2opt = p2;
               opt = quality;
               ok = true;
            }
         }
   c->p [0] = p1opt;
   c->p [1] = p2opt;
   c->p [2] = 0;
   c->s = 1;
   c->e = 1;
   return ok;
}

/*****************************************************************************/
/*                                                                           */
/* functions handling files                                                  */
/*                                                                           */
/*****************************************************************************/

void cm_class_write (cm_class_t c)
   /* writes the class polynomial to the file                                */
   /* CM_CLASS_DATADIR + "/cp_" + d + "_" + invariant + "_" + paramstr       */
   /* + ".dat"                                                               */

{
   char filename [400];
   FILE *f;
   int i;

   sprintf (filename, "%s/cp_%"PRIicl"_%c_%s.dat", CM_CLASS_DATADIR, -c.d,
            c.invariant, c.paramstr);

   if (!cm_file_open_write (&f, filename))
      exit (1);

   fprintf (f, "%"PRIicl"\n", -c.d);
   fprintf (f, "%c\n", c.invariant);
   fprintf (f, "%s\n", c.paramstr);
   fprintf (f, "%i\n", c.minpoly_deg);
   for (i = c.minpoly_deg; i >= 0; i--) {
      mpz_out_str (f, 10, c.minpoly [i]);
      if (c.field == CM_FIELD_COMPLEX) {
         fprintf (f, " ");
         mpz_out_str (f, 10, c.minpoly_complex [i]);
      }
      fprintf (f, "\n");
   }

   cm_file_close (f);
}

/*****************************************************************************/

bool cm_class_read (cm_class_t c)
   /* reads the class polynomial from a file written by cm_class_write       */
   /* If an error occurs, the return value is false.                         */

{
   char filename [400];
   FILE* f;
   int i;
   char inv;
   char pars [255];
   int_cl_t disc;

   sprintf (filename, "%s/cp_%"PRIicl"_%c_%s.dat", CM_CLASS_DATADIR, -c.d,
            c.invariant, c.paramstr);

   if (!cm_file_open_read (&f, filename))
      return false;

   if (!fscanf (f, "%"SCNicl"\n", &disc))
      return false;
   if (-disc != c.d) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** discriminant %"PRIicl" instead of %"PRIicl"\n", -disc, c.d);
      exit (1);
   }
   if (!fscanf (f, "%c", &inv))
      return false;
   if (inv != c.invariant) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** invariant '%c' instead of '%c'\n", inv, c.invariant);
      exit (1);
   }
   if (!fscanf (f, "%254s", pars))
      return false;
   if (strcmp (pars, c.paramstr)) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** parameter %s instead of %s\n", pars, c.paramstr);
      exit (1);
   }
   if (!fscanf (f, "%i", &i))
      return false;
   if (i != c.minpoly_deg) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** degree %i instead of %i\n", i, c.minpoly_deg);
      exit (1);
   }

   for (i = c.minpoly_deg; i >= 0; i--) {
      mpz_inp_str (c.minpoly [i], f, 10);
      if (c.field == CM_FIELD_COMPLEX)
         mpz_inp_str (c.minpoly_complex [i], f, 10);
   }

   cm_file_close (f);

   return true;
}

/*****************************************************************************/

static void write_conjugates (cm_class_t c, ctype *conjugate, int *conj)
   /* writes the conjugates to the file
      CM_CLASS_TMPDIR + "/tmp_" + d + "_" + invariant + "_" + paramstr + "_"
      + prec + "_conjugates.dat" */

{
   char filename [400];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%c_%s_%i_conjugates.dat",
      CM_CLASS_TMPDIR, -c.d, c.invariant, c.paramstr,
      (int) cget_prec (conjugate [0]));

   if (!cm_file_open_write (&f, filename))
      exit (1);

   for (i = 0; i < c.h; i++)
      if (conj [i] >= i) {
         cout_str (f, 16, 0, conjugate [i]);
         fprintf (f, "\n");
      }

   cm_file_close (f);
}

/*****************************************************************************/

static bool read_conjugates (cm_class_t c, ctype *conjugate, int *conj)
   /* reads the conjugates from a file written by write_conjugates
      If the file could not be openend, the return value is false. */
{
   char filename [400];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%c_%s_%i_conjugates.dat",
      CM_CLASS_TMPDIR, -c.d, c.invariant, c.paramstr,
      (int) cget_prec (conjugate [0]));

   if (!cm_file_open_read (&f, filename))
      return false;

   for (i = 0; i < c.h; i++)
      if (conj [i] >= i)
         cinp_str (conjugate [i], f, NULL, 16);

   cm_file_close (f);

   return true;
}

/*****************************************************************************/
/*                                                                           */
/* valuation at infinity and height                                          */
/*                                                                           */
/*****************************************************************************/

static double class_get_valuation (cm_class_t c)
   /* returns the (negative) order of the modular function at infinity       */

{
   double result;

   switch (c.invariant) {
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
      if (c.p [0] == 47)
         result = 1.0 / 24;
      else if (c.p [0] == 59)
         result = 1.0 / 30;
      else if (c.p [0] == 71)
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
      for (i = 0; c.p [i] != 0; i++) {
         num *= c.p [i] - 1;
         den *= c.p [i] + 1;
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
      result = (c.p [0] - 1) / (double) (24 * (c.p [0] + 1));
      break;
   default: /* should not occur */
      printf ("class_get_valuation called for unknown class ");
      printf ("invariant\n");
      exit (1);
   }

   result *= c.e;

   return result;
}

/*****************************************************************************/

static int class_get_height (cm_class_t c)
   /* in the real case, returns the binary length of the largest             */
   /* coefficient of the minimal polynomial                                  */
   /* in the complex case, returns the binary length of the largest          */
   /* coefficient with respect to the decomposition over an integral basis   */

{
   int   i, height, cand;

   height = -1;
   for (i = 0; i < c.minpoly_deg; i++) {
      cand = mpz_sizeinbase (c.minpoly [i], 2);
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

static void compute_nsystem (cm_form_t *nsystem, int *conj, cm_class_t *c,
   cm_classgroup_t cl, bool verbose)
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
         b0 = c->d % 2;
         while (!ok) {
            b0 += 2;
            if ((b0*b0 - c->d) % (4*c->p[0]) != 0)
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
      if (c->d % 2 == 0)
         neutral.b = 0;
      else
         neutral.b = 1;

      switch (c->invariant) {
         case CM_INVARIANT_J:
            b0 = c->d % 2;
            /* An even N makes c even if 2 is split so that during the
               evaluation of eta(z/2) for f1 all forms can be taken from
               the precomputed ones. */
            N = 2;
            break;
         case CM_INVARIANT_GAMMA2:
            b0 = 3 * (c->d % 2);
            /* Use an even N as for j. */
            N = 6;
            break;
         case CM_INVARIANT_GAMMA3:
            b0 = 1;
            N = 2;
            break;
         case CM_INVARIANT_ATKIN:
            N = c->p [0];
            if (c->d % 2 == 0)
               b0 = 0;
            else
               b0 = 1;
            while ((b0*b0 - c->d) % N != 0)
               b0 += 2;
            neutral.a = N;
            neutral.b = -b0;
            cm_classgroup_reduce (&neutral, c->d);
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
            if (c->d % 2 == 0)
               b0 = 2;
            else
               b0 = 1;
            while (true) {
               C = (b0*b0 - c->d) / 4;
               if (C % N == 0 && cm_nt_gcd (C / N, N) == 1)
                  break;
               b0 += 2;
            }
            neutral.a = N;
            neutral.b = -b0;
            cm_classgroup_reduce (&neutral, c->d);
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

   for (i = 0; i < cl.h; i++) {
      nsystem [i] = cl.form [i];
      conj [i] = -1;
   }

   for (i = 0; i < cl.h; i++)
      /* Pair forms yielding complex conjugate roots. */
      if (c->field == CM_FIELD_REAL) {
         if (conj [i] == -1) {
            /* form did not yet occur in a pair; look for its inverse
               with respect to neutral_class */
            nsystem [i].b = -nsystem [i].b;
            cm_classgroup_compose (&inverse, neutral, nsystem [i], c->d);
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
   for (i = 0; i < cl.h; i++)
      if (conj [i] >= i)
         correct_nsystem_entry (&(nsystem [i]), N, b0, c->d);

   /* Compute h1 and h2, only for printing. */
   if (verbose) {
      h1 = 0;
      h2 = 0;
      for (i = 0; i < cl.h; i++)
         if (conj [i] == i)
            h1++;
         else if (conj [i] > i)
            h2++;
      printf ("h = %i, h1 = %i, h2 = %i\n", c->h, h1, h2);
   }
}

/*****************************************************************************/

static fprec_t compute_precision (cm_class_t c, cm_classgroup_t cl,
   bool verbose) {
   /* returns an approximation of the required precision */

   const double C = 2114.567;
   const double pisqrtd = 3.14159265358979323846 * sqrt ((double) (-cl.d));
   const double cf = class_get_valuation (c);
   double x, binom = 1.0, simpleprec = 0, prec = 0, M;
   int_cl_t amax;
   int i, m;
   fprec_t precision;

   /* heuristic formula: log (height) = pi * sqrt (|d|) * \sum 1/A */
   for (i = 0; i < cl.h; i++)
      simpleprec += 1.0 / cl.form [i].a;
   simpleprec = ceil (simpleprec * pisqrtd / log (2.0) * cf);

   /* formula of Lemma 8 of [Sutherland11] */
   amax = 0;
   for (i = 0; i < cl.h; i++) {
      x = pisqrtd / cl.form [i].a;
      if (x < 42)
         M = log (exp (x) + C);
      else /* prevent overflow in exponential without changing the result */
         M = x;
      prec += M;
      if (cl.form [i].a > amax)
         amax = cl.form [i].a;
   }
   M = exp (pisqrtd / amax) + C;
   m = (int) ((cl.h + 1) / (M + 1));
   for (i = 1; i <= m; i++)
      binom *= (double) (cl.h - 1 + i) / i / M;
   prec = ceil ((prec + log (binom)) / log (2.0) * cf);

   if (verbose) {
      printf ("Heuristic precision bound:      %ld\n", (long int) simpleprec);
      printf ("Less heuristic precision bound: %ld\n", (long int) prec);
   }

   if (c.invariant == CM_INVARIANT_GAMMA3) {
      /* increase height estimate by bit size of sqrt (|D|)^h in constant    */
      /* coefficient                                                         */
      prec += (int) (log ((double) (-cl.d)) / log (2.0) / 2.0 * cl.h);
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

static void eval (cm_class_t c, cm_modclass_t mc, ctype rop, cm_form_t Q)

{
   switch (c.invariant) {
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
      cm_modclass_atkinhecke_level_eval_quad (mc, rop, Q.a, Q.b, c.p [0]);
      break;
   case CM_INVARIANT_SIMPLEETA:
   case CM_INVARIANT_DOUBLEETA:
   case CM_INVARIANT_MULTIETA:
         cm_modclass_multieta_eval_quad (mc, rop, Q.a, Q.b, c.p, c.e);
      break;
   case CM_INVARIANT_WEBER:
      if (c.p [0] == 1) {
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 2);
         cmul_fr (rop, rop, mc.sqrt2_over2);
      }
      else if (c.p [0] == 3)
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 1);
      else if (c.p [0] == 5) {
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 2);
         csqr (rop, rop);
         cdiv_ui (rop, rop, 2ul);
      }
      else if (c.p [0] == 7) {
         cm_modclass_f_eval_quad (mc, rop, Q.a, Q.b, 1);
         cmul_fr (rop, rop, mc.sqrt2_over2);
      }
      else if (c.p [0] == 2 || c.p [0] == 6) {
         cm_modclass_f1_eval_quad (mc, rop, Q.a, Q.b, 1);
         csqr (rop, rop);
         cmul_fr (rop, rop, mc.sqrt2_over2);
      }
      else {
         /* c.p [0] == 4 */
         cm_modclass_f1_eval_quad (mc, rop, Q.a, Q.b, 1);
         cpow_ui (rop, rop, 4ul);
         cmul_fr (rop, rop, mc.sqrt2_over4);
      }

      if (c.d % 3 == 0)
         cpow_ui (rop, rop, 3ul);

      if (c.p [0] != 3 && c.p [0] != 5)
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
   int *conj, cm_class_t c, cm_modclass_t mc, bool verbose)
   /* computes the conjugates of the singular value over Q */

{
   int i;

   for (i = 0; i < c.h; i++) {
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

static void get_int (mpz_ptr z, mpfr_ptr c)
{
   if (!cm_nt_fget_z (z, c)) {
       printf ("***** Error in rounding: ");
       fout_str (stdout, 10, 0ul, c);
       printf ("\n");
       mpz_out_str (stdout, 10, z);
       printf ("\n");
       exit (1);
   }
}

/*****************************************************************************/

static bool get_quadratic (mpz_t out1, mpz_t out2, ctype in, int_cl_t d)
   /* tries to round the complex number to an integer in the quadratic order */
   /* of discriminant d by decomposing it over the integral basis [1, omega] */
   /* with omega = sqrt (d/4) or omega = (1 + sqrt (d)) / 2                  */
   /* The return value reflects the success of the operation.                */
   /* out1 and out2 are changed.                                             */

{
   ftype   omega_i, tmp;
   bool     div4, ok;

   finit (omega_i, cget_prec (in));
   finit (tmp, cget_prec (in));

   div4 = (cm_classgroup_mod (d, (uint_cl_t) 4) == 0);
   fsqrt_ui (omega_i, -d);
   fdiv_2ui (omega_i, omega_i, 1ul);

   fdiv (tmp, in->im, omega_i);
   ok = cm_nt_fget_z (out2, tmp);

   if (ok) {
      if (div4)
         fset (tmp, in->re);
      else {
         fset_z (tmp, out2);
         fdiv_2ui (tmp, tmp, 1ul);
         fsub (tmp, in->re, tmp);
      }
      ok = cm_nt_fget_z (out1, tmp);
   }

   fclear (omega_i);
   fclear (tmp);

   return ok;
}

/*****************************************************************************/

static void round_and_print_tower (cm_class_t c, mpfrx_tower_srcptr t)
/* TODO: This function is part of a proof of concept solution; in the
   longer run, the code should be made more modular. */
{
   mpz_t z;
   int i, j, k;
   int l = t->levels;
   int *d = t->d;
   mpfrx_srcptr f;

   mpz_init (z);

   f = t->W [0][0];
   mpz_set_ui (c.tower [0][0][mpfrx_get_deg (f)], 1);
   printf ("f1 = x1^%i", mpfrx_get_deg (f));
   for (k = mpfrx_get_deg (f) - 1; k >= 0; k--) {
      get_int (z, mpfrx_get_coeff (f, k));
      mpz_set (c.tower [0][0][k], z);
      if (mpz_cmp_ui (z, 0) > 0)
         printf ("+");
      mpz_out_str (stdout, 10, z);
      if (k >= 2)
         printf ("*x1^%i", k);
      else if (k == 1)
         printf ("*x1");
   }
   printf (";\n");
   for (i = 1; i < l; i++) {
      printf ("f%i = ", i+1);
      for (j = d [i]; j >= 0; j--) {
         f = t->W [i][j];
         if (j < d [i])
            printf ("+");
         printf ("(");
         for (k = mpfrx_get_deg (f); k >= 0; k--) {
            get_int (z, mpfrx_get_coeff (f, k));
            mpz_set (c.tower [i][j][k], z);
            if (k < mpfrx_get_deg (f) && mpz_cmp_ui (z, 0) > 0)
               printf ("+");
            mpz_out_str (stdout, 10, z);
            if (k >= 2)
               printf ("*x%i^%i", i, k);
            else if (k == 1)
               printf ("*x%i", i);
         }
         if (j >= 2)
            printf (")*x%i^%i", i+1, j);
         else if (j == 1)
            printf (")*x%i", i+1);
         else
            printf (")");
      }
      printf (";\n");
   }

   mpz_clear (z);
}

/*****************************************************************************/

static void round_and_print_complex_tower (cm_class_t c, mpcx_tower_srcptr t)
/* TODO: This function is part of a proof of concept solution; in the
   longer run, the code should be made more modular. */
{
   mpz_t z1, z2;
   int i, j, k;
   int l = t->levels;
   int *d = t->d;
   mpcx_srcptr f;

   mpz_init (z1);
   mpz_init (z2);

   f = t->W [0][0];
   mpz_set_ui (c.tower [0][0][mpfrx_get_deg (f)], 1);
   mpz_set_ui (c.tower_complex [0][0][mpfrx_get_deg (f)], 0);
   printf ("f1 = x1^%i", mpcx_get_deg (f));
   for (k = mpcx_get_deg (f) - 1; k >= 0; k--) {
      get_quadratic (z1, z2, mpcx_get_coeff (f, k), c.d);
      mpz_set (c.tower [0][0][k], z1);
      mpz_set (c.tower_complex [0][0][k], z2);
      printf ("+ (");
      mpz_out_str (stdout, 10, z1);
      printf ("+(");
      mpz_out_str (stdout, 10, z2);
      printf ("*omega))");
      if (k >= 2)
         printf ("*x1^%i", k);
      else if (k == 1)
         printf ("*x1");
   }
   printf (";\n");
   for (i = 1; i < l; i++) {
      printf ("f%i = ", i+1);
      for (j = d [i]; j >= 0; j--) {
         f = t->W [i][j];
         if (j < d [i])
            printf ("+");
         printf ("(");
         for (k = mpfrx_get_deg (f); k >= 0; k--) {
            get_quadratic (z1, z2, mpcx_get_coeff (f, k), c.d);
            mpz_set (c.tower [i][j][k], z1);
            mpz_set (c.tower_complex [i][j][k], z2);
            printf ("+(");
            mpz_out_str (stdout, 10, z1);
            printf ("+(");
            mpz_out_str (stdout, 10, z2);
            printf ("*omega))");
            if (k >= 2)
               printf ("*x%i^%i", i, k);
            else if (k == 1)
               printf ("*x%i", i);
         }
         if (j >= 2)
            printf (")*x%i^%i", i+1, j);
         else if (j == 1)
            printf (")*x%i", i+1);
         else
            printf (")");
      }
      printf (";\n");
   }

   mpz_clear (z1);
   mpz_clear (z2);
}

/*****************************************************************************/

void cm_class_compute_minpoly (cm_class_t c, bool tower, bool checkpoints,
   bool disk, bool print, bool verbose)
   /* tower indicates whether the class polynomial should be decomposed as
      a Galois tower.
      checkpoints indicates whether intermediate results are to be kept in
      and potentially read from files.
      disk indicates whether the result should be written to disk.
      print indicates whether the result should be printed on screen. */
{
   cm_classgroup_t cl, cl2;
   cm_form_t *nsystem;
   int *conj;
   fprec_t prec;
   cm_modclass_t mc;
   ctype *conjugate;
   mpfrx_tower_t t;
   mpcx_tower_t tc;
   int i;
   cm_timer  clock_global, clock_local;

   cm_timer_start (clock_global);

   cm_timer_start (clock_local);
   cm_classgroup_init (&cl, c.d, verbose);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for class group: %.1f\n", cm_timer_get (clock_local));

   if (   (c.invariant == CM_INVARIANT_WEBER
           && (   (c.p [0] / 10) % 10 == 1
               || c.p [0] % 10 == 3 || c.p [0] % 10 == 4 || c.p [0] % 10 == 7))
       || (   (c.invariant == CM_INVARIANT_J || c.invariant == CM_INVARIANT_GAMMA2)
           && c.d % 4 == 0
           && ((c.d / 4) % 4 == 0 || ((c.d / 4) - 1) % 4 == 0))) {
      /* also compute class group for order of conductor less by a factor    */
      /* of 2                                                                */
      cm_timer_start (clock_local);
      cm_classgroup_init (&cl2, c.d / 4, verbose);
      cm_timer_stop (clock_local);
      if (verbose)
         printf ("--- Time for class group2: %.1f\n", cm_timer_get (clock_local));
   }
   else
      cl2.d = 0;

   nsystem = (cm_form_t *) malloc (c.h * sizeof (cm_form_t));
   conj = (int *) malloc (c.h * sizeof (int));
   cm_timer_start (clock_local);
   compute_nsystem (nsystem, conj, &c, cl, verbose);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for N-system: %.1f\n", cm_timer_get (clock_local));
   if (c.invariant == CM_INVARIANT_WEBER && c.p [0] % 100 == 17)
      prec = compute_precision (c, cl2, verbose);
   else
      prec = compute_precision (c, cl, verbose);

   conjugate = (ctype *) malloc (c.h * sizeof (ctype));
   for (i = 0; i < c.h; i++)
      if (conj [i] >= i)
         cinit (conjugate [i], prec);
   cm_timer_start (clock_local);
   if (!checkpoints || !read_conjugates (c, conjugate, conj)) {
      cm_modclass_init (&mc, cl, cl2, prec, checkpoints, verbose);
      compute_conjugates (conjugate, nsystem, conj, c, mc, verbose);
      if (checkpoints)
         write_conjugates (c, conjugate, conj);
      cm_modclass_clear (&mc);
   }
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for conjugates: %.1f\n", cm_timer_get (clock_local));

   cm_timer_start (clock_local);
   if (c.field == CM_FIELD_REAL) {
      real_compute_minpoly (c, conjugate, conj, print);
         /* Printing the full polynomial could be made optional when the
            tower is computed. */
      if (tower) {
         mpfrx_tower_init (t, c.tower_levels, c.tower_d, prec);
         mpfrcx_tower_decomposition (t, conjugate, conj);
         /* Printing should be optional here, and there should also be a
            possibility to save the field tower on disk. */
         round_and_print_tower (c, t);
      }
   }
   else {
      complex_compute_minpoly (c, conjugate, print);
      if (tower) {
         mpcx_tower_init (tc, c.tower_levels, c.tower_d, prec);
         mpcx_tower_decomposition (tc, conjugate);
         round_and_print_complex_tower (c, tc);
      }
   }
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for polynomial reconstruction: %.1f\n",
              cm_timer_get (clock_local));

   for (i = 0; i < c.h; i++)
      if (conj [i] >= i)
         cclear (conjugate [i]);
   free (conjugate);
   free (nsystem);
   free (conj);
   cm_classgroup_clear (&cl);

   cm_timer_stop (clock_global);
   if (verbose) {
      printf ("--- Total time for minimal polynomial: %.1f\n",
            cm_timer_get (clock_global));
      printf ("Height of minimal polynomial: %d\n", class_get_height (c));
   }
   if (disk)
      cm_class_write (c);
}

/*****************************************************************************/

static void real_compute_minpoly (cm_class_t c, ctype *conjugate, int *conj,
   bool print)
   /* computes the minimal polynomial of the function over Q */

{
   mpfrx_t mpol;
   fprec_t prec;
   int i;

   for (i = 0; conj [i] < i; i++);
   prec = cget_prec (conjugate [i]);
   mpfrx_init (mpol, c.minpoly_deg + 1, prec);
   mpfrcx_reconstruct_from_roots (mpol, conjugate, conj, c.minpoly_deg);

   /* the minimal polynomial is now in mpol, rounding to integral polynomial */
   for (i = 0; i < c.minpoly_deg; i++)
      if (!cm_nt_fget_z (c.minpoly [i], mpol->coeff [i])) {
         printf ("*** Accuracy not sufficient for coefficient of X^%d = ", i);
         fout_str (stdout, 10, 0ul, mpol->coeff [i]);
         printf ("\n");
         exit (1);
      }
   mpz_set_ui (c.minpoly [c.minpoly_deg], 1ul);

   if (print) {
      printf ("x^%i", c.minpoly_deg);
      for (i = c.minpoly_deg - 1; i >= 0; i--) {
         int cmp = mpz_cmp_ui (c.minpoly [i], 0);
         if (cmp != 0) {
            if (cmp > 0) {
               printf (" +");
               mpz_out_str (stdout, 10, c.minpoly [i]);
            }
            else
               mpz_out_str (stdout, 10, c.minpoly [i]);
            if (i >= 2)
               printf ("*x^%i", i);
            else if (i == 1)
               printf ("*x");
         }
      }
      printf ("\n");
   }

   mpfrx_clear (mpol);
}

/*****************************************************************************/

static void complex_compute_minpoly (cm_class_t c, ctype *conjugate,
   bool print)
   /* computes the minimal polynomial of the function over Q (sqrt D) */

{
   int_cl_t fund;
   int i;
   fprec_t prec;
   mpcx_t mpol;

   prec = cget_prec (conjugate [0]);
   mpcx_init (mpol, c.minpoly_deg + 1, prec);
   mpcx_reconstruct_from_roots (mpol, conjugate, c.minpoly_deg);

   /* the minimal polynomial is now in mpol; round to integral polynomial */
   fund = cm_classgroup_fundamental_discriminant (c.d);
   for (i = 0; i < c.h; i++) {
      if (!get_quadratic (c.minpoly [i], c.minpoly_complex [i],
         mpol->coeff[i], fund)) {
         printf ("*** accuracy not sufficient for coefficient of X^%d = ", i);
         cout_str (stdout, 10, 0ul, mpol->coeff [i]);
         printf ("\n");
         exit (1);
      }
   }

   mpcx_clear (mpol);

   if (print) {
      printf ("Minimal polynomial:\nx^%i", c.minpoly_deg);
      for (i = c.minpoly_deg - 1; i >= 0; i--) {
         printf (" + (");
         mpz_out_str (stdout, 10, c.minpoly [i]);
         if (mpz_cmp_ui (c.minpoly_complex [i], 0) >= 0)
            printf ("+");
         mpz_out_str (stdout, 10, c.minpoly_complex [i]);
         printf ("*omega) * x^%i", i);
      }
      printf ("\n");
   }
}

/*****************************************************************************/
/*                                                                           */
/* computing j-invariants                                                    */
/*                                                                           */
/*****************************************************************************/

static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose)
   /* returns a root of the minimal polynomial modulo P                      */
   /* root is changed                                                        */

{
   /* The local variables are needed only for the complex case. */
   int_cl_t fund;
   cm_timer clock;
   mpz_t *minpoly_P;
   int i;
   mpz_t omega, tmp;
   /* the second element of the integral basis of the maximal order */
   /* modulo P                                                      */

   if (c.field == CM_FIELD_REAL)
      cm_pari_oneroot (root, c.minpoly, c.minpoly_deg, P, verbose);
   else {
      mpz_init (omega);
      mpz_init (tmp);

      minpoly_P = (mpz_t*) malloc ((c.minpoly_deg + 1) * sizeof (mpz_t));
      for (i = 0; i <= c.minpoly_deg; i++)
         mpz_init (minpoly_P [i]);

      /* compute the second element of the integral basis modulo P */
      cm_timer_start (clock);
      fund = cm_classgroup_fundamental_discriminant (c.d);
      cm_nt_mpz_tonelli (omega, fund, P);
      cm_timer_stop (clock);
      if (verbose)
         printf ("--- Time for square root: %.1f\n", cm_timer_get (clock));
      if (fund % 4 != 0)
         mpz_add_ui (omega, omega, 1ul);
      mpz_set_ui (tmp, 2ul);
      mpz_invert (tmp, tmp, P);
      mpz_mul (omega, omega, tmp);
      mpz_mod (omega, omega, P);

      /* create the minimal polynomial modulo P */
      for (i = 0; i < c.minpoly_deg; i++) {
         mpz_mod (minpoly_P [i], c.minpoly_complex [i], P);
         mpz_mul (tmp, minpoly_P [i], omega);
         mpz_mod (minpoly_P [i], tmp, P);
         mpz_add (tmp, minpoly_P [i], c.minpoly [i]);
         mpz_mod (minpoly_P [i], tmp, P);
      }
      mpz_set_ui (minpoly_P [c.minpoly_deg], 1ul);

      cm_pari_oneroot (root, minpoly_P, c.minpoly_deg, P, verbose);

      mpz_clear (omega);
      mpz_clear (tmp);
      for (i = 0; i <= c.minpoly_deg; i++)
         mpz_clear (minpoly_P [i]);
      free (minpoly_P);
   }

   if (verbose) {
      printf ("Root: ");
      mpz_out_str (stdout, 10, root);
      printf ("\n");
   }
}

/*****************************************************************************/

mpz_t* cm_class_get_j_mod_P (int_cl_t d, char inv, mpz_t P, int *no,
   const char* modpoldir, bool readwrite, bool verbose)
   /* If readwrite is true, tries to read the polynomial from a file, and    */
   /* computes it and writes it to the file if it is not yet present.        */

{
   cm_class_t c;
   mpz_t *j;
   mpz_t root, d_mpz, tmp, tmp2;
   cm_timer clock;

   cm_class_init (&c, d, inv, verbose);
   if (!readwrite || !cm_class_read (c))
      cm_class_compute_minpoly (c, false, false, readwrite, false, verbose);

   cm_timer_start (clock);
   mpz_init (root);
   if (inv != CM_INVARIANT_WEBER || c.p [0] != 3 || c.p [1] != 1)
      /* avoid special case of Weber polynomial factoring over extension */
      /* of degree 3; handled in weber_cm_get_j_mod_P                    */
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
         j = weber_cm_get_j_mod_P (c, root, P, no, verbose);
         break;
      case CM_INVARIANT_DOUBLEETA:
#if 0
      case CM_INVARIANT_RAMIFIED:
#endif
         if (c.s != c.e)
            mpz_powm_ui (root, root, (unsigned long int) (c.s / c.e), P);
         j = cm_get_j_mod_P_from_modular (no, modpoldir, CM_MODPOL_DOUBLEETA,
            c.p [0] * c.p [1], root, P);
         break;
      case CM_INVARIANT_MULTIETA:
      {
         int N = 1, i;
         for (i = 0; c.p [i] != 0; i++)
            N *= c.p [i];
         if (c.s != c.e)
            mpz_powm_ui (root, root, (unsigned long int) (c.s / c.e), P);
         j = cm_get_j_mod_P_from_modular (no, modpoldir,
            CM_MODPOL_MULTIETA, N, root, P);
         break;
      }
      case CM_INVARIANT_ATKIN:
         j = cm_get_j_mod_P_from_modular (no, modpoldir, CM_MODPOL_ATKIN,
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

static mpz_t* cm_get_j_mod_P_from_modular (int *no, const char* modpoldir,
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

static mpz_t* weber_cm_get_j_mod_P (cm_class_t c, mpz_t root, mpz_t P, int *no,
   bool verbose)

{
   mpz_t* j = (mpz_t*) malloc (sizeof (mpz_t));
   mpz_t f24, tmp;
   mpfp_t one;

   mpz_init (j [0]);
   mpz_init (f24);
   mpz_init (tmp);

   if (c.p [0] != 3 || c.p [1] != 1) {
      if (c.d % 3 == 0)
         mpz_powm_ui (f24, root, 2ul, P);
      else
         mpz_powm_ui (f24, root, 6ul, P);

      if (c.p [0] != 7 || c.p [1] != 1) {
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
      }
      else {
         /* d/4 = 1 (mod 8), where d/4 is the real complex multiplication */
         /* discriminant                                                  */
         mpz_invert (tmp, f24, P);
         mpz_powm_ui (f24, tmp, 4ul, P);
         mpz_sub_ui (j [0], f24, 16ul);
         mpz_powm_ui (j [0], j [0], 3ul, P);
         mpz_invert (tmp, f24, P);
         mpz_mul (j [0], j [0], tmp);
         mpz_mod (j [0], j [0], P);
      }
   }
   else {
      /* d/4 = 5 (mod 8), we need to factor over an extension of degree 3 */
      mpz_t factor [4];
      int   i;
      mpfpx_t root_poly, f24_poly, tmp_poly, tmp2_poly;
      const unsigned long int x2 [] = {0, 0, 1};
      const unsigned long int x6 [] = {0, 0, 0, 0, 0, 0, 1};

      for (i = 0; i <= 3; i++)
         mpz_init (factor [i]);
      mpfpx_type_init (P, MPFPX_KARATSUBA);
      mpfpx_init (root_poly);
      mpfpx_init (f24_poly);
      mpfpx_init (tmp_poly);
      mpfpx_init (tmp2_poly);
      mpfp_init (one);
      mpfp_set_ui (one, 1ul);

      cm_pari_onefactor (factor, c.minpoly, c.minpoly_deg, 3, P, verbose);
      mpfpx_set_z_array (root_poly, factor, 4);
      if (verbose) {
         printf ("Root: ");
         mpfpx_out (root_poly);
      }

      if (c.d % 3 == 0)
         mpfpx_set_ui_array (f24_poly, x2, 3);
      else {
         mpfpx_set_ui_array (f24_poly, x6, 7);
         mpfpx_div_r (f24_poly, f24_poly, root_poly, one);
      }

      mpfpx_invert (f24_poly, f24_poly, root_poly, one);
      mpfpx_mul_ui (f24_poly, f24_poly, 8ul);
      mpfpx_sqr (f24_poly, f24_poly);
      mpfpx_div_r (f24_poly, f24_poly, root_poly, one);
      mpfpx_sqr (f24_poly, f24_poly);
      mpfpx_div_r (f24_poly, f24_poly, root_poly, one);
      mpfpx_invert (tmp_poly, f24_poly, root_poly, one);
      mpfpx_sub_ui (f24_poly, f24_poly, 16ul);
      mpfpx_sqr (tmp2_poly, f24_poly);
      mpfpx_div_r (tmp2_poly, tmp2_poly, root_poly, one);
      mpfpx_mul (f24_poly, tmp2_poly, f24_poly);
      mpfpx_div_r (f24_poly, f24_poly, root_poly, one);
      mpfpx_mul (f24_poly, f24_poly, tmp_poly);
      mpfpx_div_r (f24_poly, f24_poly, root_poly, one);
      mpfpx_coeff_get_z (j [0], f24_poly, 0);

      for (i = 0; i <= 3; i++)
         mpz_clear (factor [i]);
      mpfpx_clear (root_poly);
      mpfpx_clear (f24_poly);
      mpfpx_clear (tmp_poly);
      mpfpx_clear (tmp2_poly);
      mpfp_clear (one);
   }

   *no = 1;

   mpz_clear (f24);
   mpz_clear (tmp);

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
/*****************************************************************************/
