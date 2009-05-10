#include "cm_class-impl.h"

static double class_get_valuation (cm_class_t c);
   /* return some value related to heights and depending on the function     */
static int class_get_height (cm_class_t c);
   /* in the real case, returns the binary length of the largest             */
   /* coefficient of the minimal polynomial                                  */
   /* in the complex case, returns the binary length of the largest          */
   /* coefficient with respect to the decomposition over an integral basis   */

static int correct_nsystem_entry (int_cl_t *a, int_cl_t *b, int_cl_t *c,
   int_cl_t N, int_cl_t b0,
   int_cl_t neutral_class_a, int_cl_t neutral_class_b, cm_class_t cl);
   /* changes (a, b, c) to obtain an N-system and returns the embedding      */
   /* number (1 or 2) corresponding to the form                              */
static void compute_nsystem_and_embedding (int_cl_t **nsystem, int *embedding,
   cm_class_t *c, cm_classgroup_t cl, bool verbose);
   /* computes and returns an N-system and whether they lead to real or      */
   /* complex conjugates: embedding[i]==2 indicates a pair of complex        */
   /* conjugates, embedding[i]==1 a real conjugate or a single complex one.  */
static mp_prec_t compute_precision (cm_class_t c, cm_classgroup_t cl,
   bool verbose);

static void eval (cm_class_t c, cm_modclass_t mc, mpc_t rop, int_cl_t a, int_cl_t b);
static void compute_conjugates (mpc_t *conjugate, int_cl_t **nsystem,
   cm_class_t c, cm_modclass_t mc, bool verbose);
   /* computes and returns the h12 conjugates                                */
static void write_conjugates (cm_class_t c, mpc_t *conjugates);
static bool read_conjugates (cm_class_t c, mpc_t *conjugates);

static bool get_quadratic (mpz_t out1, mpz_t out2, mpc_t in, int_cl_t d);

static void get_root_mod_P (cm_class_t c, mpz_t root, mpz_t P, bool verbose);
static mpz_t* cm_get_j_mod_P_from_modular (cm_class_t c, mpz_t P, int *no,
   const char* modpoldir, bool verbose);
   /* computes the possible j-values as roots of the modular polynomial      */
   /* specified by f                                                         */

/* the remaining functions are put separately as their code is quite long    */
static void real_compute_minpoly (cm_class_t c, mpc_t *conjugate, int *embedding);
static void complex_compute_minpoly (cm_class_t c, mpc_t *conjugate);
static int doubleeta_compute_parameter (int_cl_t disc, bool verbose);
static mpz_t* weber_cm_get_j_mod_P (cm_class_t c, mpz_t P, int *no,
   bool verbose);
static mpz_t* simpleeta_cm_get_j_mod_P (cm_class_t c, mpz_t P, int *no,
   bool verbose);

/*****************************************************************************/
/*                                                                           */
/* constructor and destructor                                                */
/*                                                                           */
/*****************************************************************************/

void cm_class_init (cm_class_t *c, int_cl_t disc, char inv, bool verbose)

{
   int i;

   if (disc >= 0) {
      printf ("\n*** The discriminant must be negative.\n");
      exit (1);
   }
   else if (disc % 4 != 0 && (disc - 1) % 4 != 0) {
      printf ("\n*** %"PRIicl" is not a quadratic discriminant.\n", disc);
      exit (1);
   }
   else {
      c->d = disc;
      c->invariant = inv;
      c->p = cm_class_compute_parameter (disc, inv, verbose);
      if (!(c->p))
         exit (1);
   }
   if (verbose)
      printf ("\nDiscriminant %"PRIicl", invariant %i, parameter %i\n",
               disc, inv, c->p);

   if (inv == CM_INVARIANT_SIMPLEETA)
      c->field = CM_FIELD_COMPLEX;
   else
      c->field = CM_FIELD_REAL;

   c->h = cm_classgroup_h (NULL, NULL, disc);
#if 0
   if (inv == CM_INVARIANT_RAMIFIED)
      c->minpoly_deg = c->h / 2;
   else
#endif
      c->minpoly_deg = c->h;
   c->minpoly = (mpz_t *) malloc (c->minpoly_deg * sizeof (mpz_t));
   for (i = 0; i < c->minpoly_deg; i++)
      mpz_init (c->minpoly [i]);
   if (c->field == CM_FIELD_COMPLEX) {
      c->minpoly_complex = (mpz_t *) malloc (c->minpoly_deg * sizeof (mpz_t));
      for (i = 0; i < c->minpoly_deg; i++)
         mpz_init (c->minpoly_complex [i]);
   }
}

/*****************************************************************************/

void cm_class_clear (cm_class_t *c)

{
   int i;

   for (i = 0; i < c->minpoly_deg; i++)
      mpz_clear (c->minpoly [i]);
   free (c->minpoly);
   if (c->field == CM_FIELD_COMPLEX) {
      for (i = 0; i < c->minpoly_deg; i++)
         mpz_clear (c->minpoly_complex [i]);
      free (c->minpoly_complex);
   }
}

/*****************************************************************************/

int cm_class_compute_parameter (int_cl_t disc, int inv, bool verbose)
      /* tests whether the discriminant is suited for the chosen invariant and  */
      /* in this case computes and returns the parameter p                      */
      /* otherwise, returns 0                                                   */

{
   if (inv == CM_INVARIANT_SIMPLEETA)
   {
      /* we compute (eta (alpha/l) / eta(alpha))^e, while the modular        */
      /* equation is given for the power s                                   */
      /* We let p = 10000*l + 100*e + s                                      */
      long int l, s, e;

      l = 28;
      s = 8;
      e = 2;

      if (cm_classgroup_kronecker (disc, (int_cl_t) l) == -1) {
         if (verbose)
            printf ("*** Unsuited discriminant\n\n");
         return 0;
      }

      if (verbose)
         printf ("l = %ld, s = %ld, e = %ld\n", l, s, e);

      return (10000*l + 100*s + e);
   }
   else
   {
      switch (inv)
      {
         case CM_INVARIANT_J:
            return 1;
         case CM_INVARIANT_GAMMA2:
            if (disc % 3 == 0) {
               if (verbose) {
                  printf ("\n*** %"PRIicl" is divisible by 3, so that gamma2 ",
                          disc);
                  printf ("cannot be used.\n");
               }
               return 0;
            }
            return 1;
         case CM_INVARIANT_GAMMA3:
            if (disc % 2 == 0) {
               if (verbose) {
                  printf ("\n*** %"PRIicl" is divisible by 4, so that gamma3 ",
                          disc);
                  printf ("cannot be used.\n");
               }
               return 0;
            }
            return 1;
         case CM_INVARIANT_WEBER:
            if (disc % 32 == 0) {
               if (verbose) {
                  printf ("\n*** %"PRIicl" is divisible by 32, so that the Weber ",
                          disc);
                  printf ("functions cannot be used.\n");
               }
               return 0;
            }

            /* p is what I used to call m; if disc is 1 or 5 mod 8, I add 10 */
            if ((disc - 5) % 8 == 0)
               return 15;
            else if ((disc - 1) % 8 == 0)
               return 11;
            else if ((disc/4) % 8 < 0)
               return ((disc/4) % 8 + 8);
            else
               return ((disc/4) % 8);
         case CM_INVARIANT_DOUBLEETA:
            return doubleeta_compute_parameter (disc, verbose);
#if 0
         case CM_INVARIANT_RAMIFIED:
            return ramified_compute_parameter (disc, verbose);
#endif
         default: /* should not occur */
            printf ("class_compute_parameter called for "
                    "unknown class invariant %d\n", inv);
            exit (1);
      }

      return 0;
   }
}

/*****************************************************************************/

static int doubleeta_compute_parameter (int_cl_t disc, bool verbose)
   /* we compute                                                             */
   /* eta (alpha/p1)*eta(alpha/p2) / eta(alpha)*eta(alpha/(p1*p2))           */
   /* p1, p2 >= 3 are determined such that they are splitting or ramified,   */
   /* not dividing the conductor and 24 | (p1-1)(p1-1), i. e. either         */
   /* p1 \equiv 1 mod 3 and p2 \equiv 1 mod 4,                               */
   /* or p1 \equiv 1 mod 12 and p2 arbitrary.                                */
   /* We let p2 >= p1 and put p = 100*p2 + p1.                               */
   /* The return value is 0 when the product p1*p2 is larger than 2000,      */
   /* for which the modular polynomials have not been computed yet.          */

{
   unsigned long int prime = 2;
   int mod1 = 0, mod3 = 0, mod4 = 0, mod12 = 0, mod1_2 = 0;
   /* the smallest suitable primes congruent to 1 mod 1, 3, 4 or 12;      */
   /* 0 if no such prime was found so far. The strategy is to choose      */
   /* the better of mod1*mod12 and mod3*mod4 for p1*p2.                   */
   /* mod1_2 is the second smallest such prime. It is needed if the first */
   /* non-inert prime is 1 modulo 12 and divides the discriminant. Then   */
   /* mod1 = mod12 = mod3 = mod4, and we would choose this ramified prime */
   /* twice. Thus, we have to use mod12*mod1_2.                           */
   int p1, p2, tmp;

   while (mod1 == 0 || mod3 == 0 || mod4 == 0 || mod12 == 0 || mod1_2 == 0)
   {
      prime = cm_nt_next_prime (prime);
      while (cm_classgroup_kronecker (disc, (int_cl_t) prime) == -1
             || disc % ((int_cl_t) (prime*prime)) == 0)
         prime = cm_nt_next_prime (prime);

      if (mod1 == 0)
         mod1 = prime;
      else if (mod1_2 == 0)
         mod1_2 = prime;
      if ((prime - 1) % 3 == 0 && mod3 == 0)
         mod3 = prime;
      if ((prime - 1) % 4 == 0 && mod4 == 0)
         mod4 = prime;
      if ((prime - 1) % 12 == 0 && mod12 == 0)
         mod12 = prime;
   }

   if (mod12 == mod1 && cm_classgroup_kronecker (disc, (int_cl_t) mod1) == 0) {
      if (verbose) {
         printf ("Possible value for (p1,p2) for d = %"PRIicl": ", disc);
         printf ("(%d, %d)\n", mod1_2, mod12);
      }
      p1 = mod1_2;
      p2 = mod12;
   }
   else {
      if (verbose) {
         printf ("Possible values for (p1,p2) for d = %"PRIicl": ", disc);
         printf ("(%d, %d) and (%d, %d)\n", mod1, mod12, mod3, mod4);
      }
      if ( (mod1-1) * (mod12-1) / (double) ((mod1+1) * (mod12+1))
            <(mod3-1) * (mod4-1) / (double) ((mod3+1) * (mod4+1)))
      {
         p1 = mod1;
         p2 = mod12;
      }
      else
      {
         p1 = mod3;
         p2 = mod4;
      }
   }

   if (p1 > p2)
   {
      tmp = p1;
      p1 = p2;
      p2 = tmp;
   }

   if (verbose) {
      printf ("p1 = %d\n", p1);
      printf ("p2 = %d\n", p2);
   }

   if (p1*p2 > 2000)
   {
      if (verbose) {
         printf ("\n*** The modular polynomial for this parameter combination ");
         printf ("does not exist yet.\n");
      }
      return 0;
   }

   return 100*p2 + p1;
}

/*****************************************************************************/
/*                                                                           */
/* functions handling files                                                  */
/*                                                                           */
/*****************************************************************************/

void cm_class_write (cm_class_t c)
   /* writes the class polynomial to the file                                */
   /* CLASS_DATADIR + "/cp_" + d + "_" + invariant + p + ".dat"              */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/cp_%"PRIicl"_%c%i.dat", CM_CLASS_DATADIR, -c.d,
            c.invariant, c.p);

   if (!cm_file_open_write (&f, filename))
      exit (1);

   fprintf (f, "%"PRIicl"\n", -c.d);
   fprintf (f, "%c\n", c.invariant);
   fprintf (f, "%i\n", c.p);
   fprintf (f, "%i\n", c.minpoly_deg);
   if (c.field == CM_FIELD_REAL)
      fprintf (f, "1\n");
   else
      fprintf (f, "1 0\n");
   for (i = c.minpoly_deg - 1; i >= 0; i--) {
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
   /* reads the class polynomial from the file                               */
   /* CM_CLASS_DATADIR + "/cp_" + d + "_" + invariant + p + ".dat"           */
   /* If an error occurs, the return value is false.                         */

{
   char filename [255];
   FILE* f;
   int i;
   char inv;
   int_cl_t disc;

   sprintf (filename, "%s/cp_%"PRIicl"_%c%i.dat", CM_CLASS_DATADIR, -c.d,
      c.invariant, c.p);

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
   if (!fscanf (f, "%i", &i))
      return false;
   if (i != c.p) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** parameter %i instead of %i\n", i, c.p);
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

   /* skip the leading 1 and possibly 0*/
   if (!fscanf (f, "%i", &i))
      return false;
   if (c.field == CM_FIELD_COMPLEX)
      if (!fscanf (f, "%i", &i))
         return false;

   for (i = c.minpoly_deg - 1; i >= 0; i--) {
      mpz_inp_str (c.minpoly [i], f, 10);
      if (c.field == CM_FIELD_COMPLEX)
         mpz_inp_str (c.minpoly_complex [i], f, 10);
   }

   cm_file_close (f);

   return true;
}

/*****************************************************************************/

static void write_conjugates (cm_class_t c, mpc_t *conjugate)
   /* writes the conjugates to the file                                      */
   /* CM_CLASS_TMPDIR + "/tmp_" + d + "_" + invariant + p + "_" + prec +     */
   /* "_conjugates.dat"                                                      */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%c%i_%i_conjugates.dat", CM_CLASS_TMPDIR,
      -c.d, c.invariant, c.p, (int) mpc_get_prec (conjugate [0]));

   if (!cm_file_open_write (&f, filename))
      exit (1);

   for (i = 0; i < c.h12; i++) {
      mpc_out_str (f, 16, 0, conjugate [i], MPC_RNDNN);
      fprintf (f, "\n");
   }

   cm_file_close (f);
}

/*****************************************************************************/

static bool read_conjugates (cm_class_t c, mpc_t *conjugate)
   /* reads the conjugates from the file                                     */
   /* CM_CLASS_TMPDIR + "/tmp_" + d + "_" + invariant + p + "_" + prec +     */
   /* "_conjugates.dat"                                                      */
   /* If the file could not be openend, the return value is false.           */
{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%c%i_%i_conjugates.dat", CM_CLASS_TMPDIR,
      -c.d, c.invariant, c.p, (int) mpc_get_prec (conjugate [0]));

   if (!cm_file_open_read (&f, filename))
      return false;

   for (i = 0; i < c.h12; i++) {
      mpc_inp_str (conjugate [i], f, NULL, 16, MPC_RNDNN);
   }

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
   int p1, p2;
   int l, e;

   switch (c.invariant)
   {
   case CM_INVARIANT_J:
      result = 1;
      break;
   case CM_INVARIANT_GAMMA2:
      result = 1 / (double) 3;
   case CM_INVARIANT_GAMMA3:
      result = 0.5;
      break;
   case CM_INVARIANT_WEBER:
      result = 1 / (double) 72;
      if (c.p == 2 || c.p == 6 || c.p == 7)
         result *= 2;
      else if (c.p == 3 || c.p == 4)
         result *= 4;
      if (c.d % 3 == 0)
         result *= 3;
      break;
   case CM_INVARIANT_DOUBLEETA:
      p1 = c.p % 100;
      p2 = c.p / 100;
      result = ((p1-1)*(p2-1) / (double) (12 * (p1+1) * (p2+1)));
      break;
#if 0
   case CM_INVARIANT_RAMIFIED:
      /* first guess: half the precision as CM_INVARIANT_DOUBLEETA */
      p1 = c.p % 100;
      p2 = c.p / 100;
      result = ((p1-1)*(p2-1) / (double) (24 * (p1+1) * (p2+1)));
      break;
#endif
   case CM_INVARIANT_SIMPLEETA:
      l = c.p / 10000;
      e = (c.p / 100) % 100;
      return (e * (l-1) / (double) (24 * (l+1)));
      break;
   default: /* should not occur */
      printf ("class_get_valuation called for unknown class ");
      printf ("invariant\n");
      exit (1);
   }

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

static int correct_nsystem_entry (int_cl_t *a, int_cl_t *b, int_cl_t *c,
   int_cl_t N, int_cl_t b0,
   int_cl_t neutral_class_a, int_cl_t neutral_class_b, cm_class_t cl)
   /* changes the form by a unimodular transformation so that b is           */
   /* congruent to b0 modulo 2*N                                             */
   /* Furthermore, determines and returns via type how many conjugates (real */
   /* or complex conjugate ones) correspond to the form.                     */
   /* In the real case, forms whose product is equivalent to neutral_class   */
   /* correspond to conjugate complex values. The lexicographically smaller  */
   /* of two such forms gets the return value 2 to indicate that it          */
   /* corresponds to two conjugates; the other form gets return value 0,     */
   /* is in fact not corrected to the N-system condition and may be dropped. */
   /* If the square of a form equals the neutral class, then the conjugate   */
   /* is real and the return value becomes 1.                                */
   /* In the complex case, the return value is always 1.                     */
   /* For the double eta quotients in the ramified case, the rule is more    */
   /* complex and explained in detail in the code.                           */

{
   int_cl_t inverse_a, inverse_b, mu, tmp;
   int type;

#if 0
   if (cl.invariant == CM_INVARIANT_RAMIFIED)
   {
   }
   else
#endif
   if (cl.field == CM_FIELD_REAL)
   {
      /* check for the inverse of the form with respect  to */
      /* neutral_class                                      */
      cm_classgroup_compose (&inverse_a, &inverse_b,
         neutral_class_a, neutral_class_b, *a, -(*b), cl.d);
      if (*a == inverse_a && *b == inverse_b)
      /* the conjugate is real */
         type = 1;
      /* the conjugate is complex, test whether its form is the */
      /* lexicographically smaller one                          */
      else if (*a < inverse_a || (*a == inverse_a && *b < inverse_b))
         type = 2;
      else
         type = 0;
   }
   else
      type = 1;

   if (type != 0)
   {
      /* modify a, b and c by unimodular transformations    */
      /* such that gcd (a, bN) = 1 and b \equiv b0 mod 2*bN */
      while (cm_classgroup_gcd (*a, N) != 1)
      {
         mu = 0;
         tmp = *c;
         while (cm_classgroup_gcd (tmp, N) != 1 && mu < 10)
         {
            mu++;
            tmp = *c + mu*(mu*(*a) - (*b));
         }
         if (mu == 10)
         {
            mu = 1;
            tmp = *c + mu*(mu*(*a) - (*b));
         }
         *b = 2*mu*(*a) - (*b);
         *c = *a;
         *a = tmp;
      }
      mu = 0;
      while ((*b - 2*mu*(*a) - b0) % (2*N) != 0)
         mu++;
      *c += mu*(mu*(*a) - (*b));
      *b -= 2*mu*(*a);
   }

   return type;
}

/*****************************************************************************/

static void compute_nsystem_and_embedding (int_cl_t **nsystem, int *embedding,
   cm_class_t *c, cm_classgroup_t cl, bool verbose)
   /* computes an N-system, or to be more precise, some part of an N-system  */
   /* that yields all different conjugates up to complex conjugation, as     */
   /* well as the nature (real/complex) of the corresponding conjugates      */

{
   int_cl_t b0, N;
   int_cl_t neutral_class_a, neutral_class_b;
   int i;
//   int rem;

   c->h12 = 0;
   c->h1 = 0;
   c->h2 = 0;

   if (c->invariant == CM_INVARIANT_SIMPLEETA) {
      bool ok = false;
      int l, s, e;

      l = c->p / 10000;
      s = (c->p / 100) % 100;
      e = c->p % 100;
      N = l * s / e;

      if (l != 4) {
         b0 = c->d % 2;
         while (!ok) {
            b0 += 2;
            if ((b0*b0 - c->d) % (4*l) != 0)
               ok = false;
            else if (l != 2 && (s/e) % 2 == 0 && (b0 - 1) % 4 != 0)
               ok = false;
            else if (l != 2 && (s/e) % 3 == 0 && b0 % 3 != 0)
               ok = false;
            else
               ok = true;
         }
      }
      else {
#if 0
         rem  = (int) cm_classgroup_mod (c->d, (uint_cl_t) 16ul);
         if (rem == 4)
            b0 = (int_cl_t) 2;
         else if (rem == 1)
            b0 = (int_cl_t) 1;
         else /* rem == 9 */
            b0 = (int_cl_t) 3;
#endif
        b0 = (int_cl_t) -7;
      }
      if (verbose)
         printf ("l %i\ns %i\ne %i\nN %"PRIicl"\nb0 %"PRIicl"\n", l, s, e, N, b0);
   }

   else {
      neutral_class_a = 1;
      if (c->d % 2 == 0)
         neutral_class_b = 0;
      else
         neutral_class_b = 1;

      switch (c->invariant) {
         case CM_INVARIANT_J:
            b0 = c->d % 2;
            N = 1;
            break;
         case CM_INVARIANT_GAMMA2:
            b0 = 3 * (c->d % 2);
            N = 3;
            break;
         case CM_INVARIANT_GAMMA3:
            b0 = 1;
            N = 2;
            break;
         case CM_INVARIANT_WEBER:
            if (c->d % 4 != 0)
               c->d <<= 2;
            neutral_class_a = 1;
            neutral_class_b = 0;
            b0 = 0;
            N = 48;
            break;
         case CM_INVARIANT_DOUBLEETA:
#if 0
         case CM_INVARIANT_RAMIFIED:
#endif
            N = (c->p/100)*(c->p%100);
            if (c->d % 2 == 0)
               b0 = 2;
            else
               b0 = 1;
            while ((b0*b0 - c->d) % N != 0)
               b0 += 2;
            neutral_class_a = N;
            neutral_class_b = -b0;
            cm_classgroup_reduce (&neutral_class_a, &neutral_class_b, c->d);
            break;
         default: /* should not occur */
            printf ("compute_nsystem_and_embedding called for ");
            printf ("unknown class invariant\n");
            exit (1);
      }
   }

   for (i = 0; i < cl.h12; i++) {
      nsystem [c->h12][0] = cl.form [i][0];
      nsystem [c->h12][1] = cl.form [i][1];
      nsystem [c->h12][2] = cl.form [i][2];
      embedding [c->h12] = correct_nsystem_entry (&(nsystem [c->h12][0]),
         &(nsystem [c->h12][1]), &(nsystem [c->h12][2]),
         N, b0, neutral_class_a, neutral_class_b, *c);
#if 0
      if (embedding [c->h12] != 0)
         printf ("[%"PRIicl" %"PRIicl" %"PRIicl"]: %i\n", nsystem [c->h12][0],
            nsystem [c->h12][1], nsystem [c->h12][2], embedding [c->h12]);
#endif
      if (embedding [c->h12] == 1) {
         c->h1++;
         c->h12++;
      }
      else if (embedding [c->h12] == 2) {
         c->h2++;
         c->h12++;
      }
      /* possibly include the inverse form */
      if (cl.form [i][1] > 0 && cl.form [i][1] < cl.form [i][0]
          && cl.form [i][0] < cl.form [i][2]) {
         nsystem [c->h12][0] = cl.form [i][0];
         nsystem [c->h12][1] = -cl.form [i][1];
         nsystem [c->h12][2] = cl.form [i][2];
         embedding [c->h12] = correct_nsystem_entry (&(nsystem [c->h12][0]),
            &(nsystem [c->h12][1]), &(nsystem [c->h12][2]),
            N, b0, neutral_class_a, neutral_class_b, *c);
#if 0
         if (embedding [c->h12] != 0)
            printf ("[%"PRIicl" %"PRIicl" %"PRIicl"]: %i\n",
               nsystem [c->h12][0], nsystem [c->h12][1],
               nsystem [c->h12][2], embedding [c->h12]);
#endif
         if (embedding [c->h12] == 1) {
            c->h1++;
            c->h12++;
         }
         else if (embedding [c->h12] == 2) {
            c->h2++;
            c->h12++;
         }
      }
   }
   if (verbose)
      printf ("h = %i, h1 = %i, h2 = %i\n", c->h, c->h1, c->h2);
}

/*****************************************************************************/

static mp_prec_t compute_precision (cm_class_t c, cm_classgroup_t cl,
   bool verbose) {
   /* returns an approximation of the required precision */

   double precision = 0;
   int i;

   for (i = 0; i < cl.h12; i++)
      if (   cl.form [i][1] == 0 || cl.form [i][0] == cl.form [i][1]
          || cl.form [i][0] == cl.form [i][2])
         precision += 1.0 / cl.form [i][0];
      else
         precision += 2.0 / cl.form [i][0];

   /* the constant is pi / log (2)    */
   precision *= 4.5323601418
         * sqrt ((double) (-c.d)) * class_get_valuation (c);
   if (verbose)
      printf ("Estimated precision: %.1f\n", precision);
   if (c.invariant == CM_INVARIANT_WEBER)
      precision += (precision > 100 ? precision / 3 : 30);
   else
      precision += (precision > 2000 ? precision / 100 : 20);
   if (verbose)
      printf ("Final precision: %d\n", (int) precision);

   return (mp_prec_t) precision;
}

/*****************************************************************************/

static void eval (cm_class_t c, cm_modclass_t mc, mpc_t rop, int_cl_t a, int_cl_t b)

{
   mpc_t tmp, tmp2;
   int   p1, p2;

   switch (c.invariant) {
   case CM_INVARIANT_J:
      cm_modclass_j_eval_quad (mc, rop, a, b);
      break;
   case CM_INVARIANT_GAMMA2:
      cm_modclass_gamma2_eval_quad (mc, rop, a, b);
      break;
   case CM_INVARIANT_GAMMA3:
      cm_modclass_gamma3_eval_quad (mc, rop, a, b);
      break;
   case CM_INVARIANT_DOUBLEETA:
#if 0
   case CM_INVARIANT_RAMIFIED:
#endif
      p1 = c.p % 100;
      p2 = c.p / 100;

      mpc_init2 (tmp, mc.m.prec);
      mpc_init2 (tmp2, mc.m.prec);

      if (p1 == p2) {
         cm_modclass_eta_eval_quad (mc, rop, a * p1, b);
         mpc_sqr (rop, rop, MPC_RNDNN);
      }
      else {
         cm_modclass_eta_eval_quad (mc, rop, a * p1, b);
         cm_modclass_eta_eval_quad (mc, tmp, a * p2, b);
         mpc_mul (rop, rop, tmp, MPC_RNDNN);
      }
      cm_modclass_eta_eval_quad (mc, tmp, a, b);
      cm_modclass_eta_eval_quad (mc, tmp2, a * p1 * p2, b);
      mpc_mul (tmp, tmp, tmp2, MPC_RNDNN);
      mpc_div (rop, rop, tmp, MPC_RNDNN);

      mpc_clear (tmp);
      mpc_clear (tmp2);
      break;
   case CM_INVARIANT_SIMPLEETA:
      mpc_init2 (tmp, mc.m.prec);

      cm_modclass_eta_eval_quad (mc, rop, a * (c.p / 10000), b);
      cm_modclass_eta_eval_quad (mc, tmp, a, b);
      mpc_div (rop, rop, tmp, MPC_RNDNN);
      mpc_pow_ui (rop, rop, (unsigned long int) (c.p % 100));

      mpc_clear (tmp);
      break;
   case CM_INVARIANT_WEBER:
      if (c.p == 1 || c.p == 11) {
         cm_modclass_f_eval_quad (mc, rop, a, b);
         mpc_mul_fr (rop, rop, mc.sqrt2_over2, MPC_RNDNN);
      }
      else if (c.p == 2 || c.p == 6) {
         cm_modclass_f1_eval_quad (mc, rop, a, b);
         mpc_sqr (rop, rop, MPC_RNDNN);
         mpc_mul_fr (rop, rop, mc.sqrt2_over2, MPC_RNDNN);
      }
      else if (c.p == 3) {
         cm_modclass_f_eval_quad (mc, rop, a, b);
         mpc_pow_ui (rop, rop, 4ul);
         mpc_div_ui (rop, rop, 2ul, MPC_RNDNN);
      }
      else if (c.p == 4) {
         cm_modclass_f1_eval_quad (mc, rop, a, b);
         mpc_pow_ui (rop, rop, 4ul);
         mpc_mul_fr (rop, rop, mc.sqrt2_over4, MPC_RNDNN);
      }
      else if (c.p == 5 || c.p == 15)
         cm_modclass_f_eval_quad (mc, rop, a, b);
      else {
         cm_modclass_f_eval_quad (mc, rop, a, b);
         mpc_sqr (rop, rop, MPC_RNDNN);
         mpc_mul_fr (rop, rop, mc.sqrt2_over2, MPC_RNDNN);
      }

      if (c.d % 3 == 0)
         mpc_pow_ui (rop, rop, 3ul);

      if (c.p != 3 && c.p != 5 && c.p != 15)
         if (cm_classgroup_kronecker ((int_cl_t) 2, a) == -1)
            mpc_neg (rop, rop, MPC_RNDNN);

      break;
   default: /* should not occur */
      printf ("class_eval called for unknown class invariant\n");
      exit (1);
   }
}

/*****************************************************************************/

static void compute_conjugates (mpc_t *conjugate, int_cl_t **nsystem,
   cm_class_t c, cm_modclass_t mc, bool verbose)
   /* computes the conjugates of the singular value over Q */

{
   int i;

   for (i = 0; i < c.h12; i++) {
      eval (c, mc, conjugate [i], nsystem [i][0], nsystem [i][1]);
      if (verbose && i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   if (verbose)
      printf ("\n");
}

/*****************************************************************************/

void cm_class_compute_minpoly (cm_class_t c, bool checkpoints, bool write,
   bool verbose)
   /* checkpoints indicates whether intermediate results are to be kept in   */
   /* and potentially read from files.                                       */
   /* write indicates whether the result should be written to disk.          */
{
   cm_classgroup_t cl;
   int_cl_t **nsystem;
   int *embedding;
   mp_prec_t prec;
   cm_modclass_t mc;
   mpc_t *conjugate;
   cm_timer  clock_global, clock_local;
   int i;

   cm_timer_start (clock_global);

   cm_timer_start (clock_local);
   cm_classgroup_init (&cl, c.d, checkpoints, verbose);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for class group: %.1f\n", cm_timer_get (clock_local));

   nsystem = (int_cl_t **) malloc (c.h * sizeof (int_cl_t *));
   for (i = 0; i < c.h; i++)
      nsystem [i] = (int_cl_t *) malloc (3 * sizeof (int_cl_t));
   embedding = (int*) malloc (c.h * sizeof (int));
   cm_timer_start (clock_local);
   compute_nsystem_and_embedding (nsystem, embedding, &c, cl, verbose);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for N-system: %.1f\n", cm_timer_get (clock_local));
   nsystem = (int_cl_t **) realloc (nsystem, c.h12 * sizeof (int_cl_t *));
   embedding = (int *) realloc (embedding, c.h12 * sizeof (int));
   prec = compute_precision (c, cl, verbose);

   conjugate = (mpc_t *) malloc (c.h12 * sizeof (mpc_t));
   for (int i = 0; i < c.h12; i++)
      mpc_init2 (conjugate [i], prec);
   cm_timer_start (clock_local);
   if (!checkpoints || !read_conjugates (c, conjugate)) {
      cm_modclass_init (&mc, cl, prec, checkpoints, verbose);
      compute_conjugates (conjugate, nsystem, c, mc, verbose);
      if (checkpoints)
         write_conjugates (c, conjugate);
      cm_modclass_clear (&mc);
   }
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for conjugates: %.1f\n", cm_timer_get (clock_local));

   for (i = 0; i < c.h12; i++)
      free (nsystem [i]);
   free (nsystem);

   cm_timer_start (clock_local);
   if (c.field == CM_FIELD_REAL)
      real_compute_minpoly (c, conjugate, embedding);
   else
      complex_compute_minpoly (c, conjugate);
   cm_timer_stop (clock_local);
   if (verbose)
      printf ("--- Time for polynomial reconstruction: %.1f\n",
              cm_timer_get (clock_local));

   /* is done in real_compute_minpoly or complex_compute_minpoly */
   /*
   for (i = 0; i < c.h12; i++)
      mpc_clear (conjugate [i]);
   free (conjugate);
   */
   free (embedding);

   cm_timer_stop (clock_global);
   if (verbose) {
      printf ("--- Total time for minimal polynomial: %.1f\n",
            cm_timer_get (clock_global));
      printf ("Height of minimal polynomial: %d\n", class_get_height (c));
   }
   if (write)
      cm_class_write (c);
}

/*****************************************************************************/

static void real_compute_minpoly (cm_class_t c, mpc_t *conjugate, int *embedding)
   /* computes the minimal polynomial of the function over Q                 */
   /* frees conjugates                                                       */

{
   mpfrx_t *factors;
      /* holds the small factors which must be multiplied together to obtain */
      /* the minimal polynomial.                                             */
   mpfrx_t mpol;
   int left = 0, right = c.h12 - 1;
      /* the next free places in the array of factors, from the left and     */
      /* from the right                                                      */
   int i;

   /* Copy the real conjugates and the pairs of complex ones into            */
   /* polynomials over the reals.                                            */
   /* Put the real ones to the right, the complex ones to the left.          */
   /* To save memory, free the conjugates at the same time.                  */
   factors = (mpfrx_t*) malloc (c.h12 * sizeof (mpfrx_t));
   for (i = 0; i < c.h12; i++) {
      if (embedding [i] == 1) {
         mpfrx_init (factors [right], 2, mpfr_get_prec (conjugate [i]->re));
         factors [right]->deg = 1;
         mpfr_set_ui (factors [right]->coeff [1], 1ul, GMP_RNDN);
         mpfr_neg (factors [right]->coeff [0], conjugate [i]->re, GMP_RNDN);
         right--;
      }
      else {
         mpfrx_init (factors [left], 3, mpfr_get_prec (conjugate [i]->re));
         factors [left]->deg = 2;
         mpfr_set_ui (factors [left]->coeff [2], 1ul, GMP_RNDN);
         mpfr_mul_2ui (factors [left]->coeff [1], conjugate [i]->re, 1ul,
            GMP_RNDN);
         mpfr_neg (factors [left]->coeff [1], factors [left]->coeff [1],
            GMP_RNDN);
         mpc_norm (factors [left]->coeff [0], conjugate [i], GMP_RNDN);
         left++;
      }
      mpc_clear (conjugate [i]);
   }
   free (conjugate);

   /* Group the linear polynomials by pairs. */
   for (i = 0; i < c.h1/2; i++) {
      mpfrx_mul (factors [c.h2+i], factors [c.h2+i], factors [c.h12-1-i]);
      mpfrx_clear (factors [c.h12-1-i]);
   }

   mpfrx_init (mpol, c.minpoly_deg + 1, factors [0]->prec);
   mpfrx_reconstruct (mpol, factors, (c.h1+1)/2 + c.h2);
   for (i = 0; i < (c.h1+1)/2 + c.h2; i++)
      mpfrx_clear (factors [i]);
   free (factors);

   /* the minimal polynomial is now in mpol, rounding to integral polynomial */
   for (i = 0; i < c.minpoly_deg; i++)
      if (!cm_nt_mpfr_get_z (c.minpoly [i], mpol->coeff [i])) {
         printf ("*** accuracy not sufficient for coefficient of X^%d = ", i);
         mpfr_out_str (stdout, 10, 0ul, mpol->coeff [i], GMP_RNDN);
         printf ("\n");
         exit (1);
      }

/*
      printf ("x^%i", c.minpoly_deg);
      for (i = c.minpoly_deg - 1; i >= 0; i--) {
         printf (" + (");
         mpz_out_str (stdout, 10, c.minpoly [i]);
         printf (") * x^%i", i);
      }
      printf ("\n");
*/

   mpfrx_clear (mpol);
}

/*****************************************************************************/

static bool get_quadratic (mpz_t out1, mpz_t out2, mpc_t in, int_cl_t d)
   /* tries to round the complex number to an integer in the quadratic order */
   /* of discriminant d by decomposing it over the integral basis [1, omega] */
   /* with omega = sqrt (d/4) or omega = (1 + sqrt (d)) / 2                  */
   /* The return value reflects the success of the operation.                */
   /* out1 and out2 are changed.                                             */

{
   mpfr_t   omega_i, tmp;
   bool     div4, ok;

   mpfr_init2 (omega_i, mpc_get_prec (in));
   mpfr_init2 (tmp, mpc_get_prec (in));

//   printf ("Fund %"PRIicl"\n", d);
   div4 = (cm_classgroup_mod (d, (uint_cl_t) 4) == 0);
   mpfr_sqrt_ui (omega_i, -d, GMP_RNDN);
   mpfr_div_2ui (omega_i, omega_i, 1ul, GMP_RNDN);

   mpfr_div (tmp, in->im, omega_i, GMP_RNDN);
   ok = cm_nt_mpfr_get_z (out2, tmp);

   if (!ok)
      return false;
   else {
      if (div4)
         mpfr_set (tmp, in->re, GMP_RNDN);
      else {
         mpfr_set_z (tmp, out2, GMP_RNDN);
         mpfr_div_2ui (tmp, tmp, 1ul, GMP_RNDN);
         mpfr_sub (tmp, in->re, tmp, GMP_RNDN);
      }
      ok = cm_nt_mpfr_get_z (out1, tmp);
      return ok;
   }

   mpfr_clear (omega_i);
   mpfr_clear (tmp);
}

/*****************************************************************************/

static void complex_compute_minpoly (cm_class_t c, mpc_t *conjugate)
   /* computes the minimal polynomial of the function over Q                 */
   /* frees conjugates                                                       */

{
   int_cl_t fund;
   int i;
   mpcx_t *factors, mpol;
      /* holds the small factors which must be multiplied together to obtain */
      /* the minimal polynomial.                                             */

   /* copy the conjugates into polynomials over the complex numbers */
   /* To save memory, free the conjugates at the same time.         */
   factors = (mpcx_t*) malloc (c.h12 * sizeof (mpcx_t));
   for (i = 0; i < c.h12; i++) {
      mpcx_init (factors [i], 2, mpfr_get_prec (conjugate [i]->re));
      factors [i]->deg = 1;
      mpc_set_ui (factors [i]->coeff [1], 1ul, MPC_RNDNN);
      mpc_set (factors [i]->coeff [0], conjugate [i], MPC_RNDNN);
      mpc_neg (factors [i]->coeff [0], factors [i]->coeff[0], MPC_RNDNN);
      mpc_clear (conjugate [i]);
   }
   free (conjugate);

   mpcx_init (mpol, c.minpoly_deg + 1, factors [0]->prec);
   mpcx_reconstruct (mpol, factors, c.h12);
   for (i = 0; i < c.h12; i++)
      mpcx_clear (factors [i]);
   free (factors);

   /* the minimal polynomial is now in mpol */
   /* rounding to integral polynomial       */
   fund = cm_classgroup_fundamental_discriminant (c.d);
   for (i = 0; i < c.h12; i++) {
      if (!get_quadratic (c.minpoly [i], c.minpoly_complex [i],
         mpol->coeff[i], fund)) {
         printf ("*** accuracy not sufficient for coefficient of X^%d = ", i);
         mpc_out_str (stdout, 10, 0ul, mpol->coeff [i], MPC_RNDNN);
         printf ("\n");
         exit (1);
      }
   }

   mpcx_clear (mpol);

   printf ("Minimal polynomial real part:\n1");
   for (i = c.minpoly_deg - 1; i >= 0; i--) {
      printf (" ");
      mpz_out_str (stdout, 10, c.minpoly [i]);
   }
   printf ("\n");
   printf ("Minimal polynomial part in the second element of the canonical integral basis:\n0");
   for (i = c.minpoly_deg - 1; i >= 0; i--) {
      printf (" ");
      mpz_out_str (stdout, 10, c.minpoly_complex [i]);
   }
   printf ("\n");
}

/*****************************************************************************/
/*                                                                           */
/* computing j-invariants                                                    */
/*                                                                           */
/*****************************************************************************/

mpz_t* cm_class_get_j_mod_P (int_cl_t d, char inv, mpz_t P, int *no,
   const char* modpoldir, bool read, bool verbose)
   /* If read is true, tries to read the polynomial from a file. */

{
   cm_class_t c;
   mpz_t *j;
   mpz_t root, d_mpz, tmp, tmp2;
   cm_timer clock;
   int i;

   cm_class_init (&c, d, inv, verbose);
   if (!read || !cm_class_read (c))
      cm_class_compute_minpoly (c, false, false, verbose);
   cm_timer_start (clock);
   switch (inv)
   {
      case CM_INVARIANT_J:
         j = (mpz_t*) malloc (sizeof (mpz_t));
         mpz_init (j [0]);
         get_root_mod_P (c, j [0], P, verbose);
         *no = 1;
         break;
      case CM_INVARIANT_GAMMA2:
         j = (mpz_t*) malloc (sizeof (mpz_t));
         mpz_init (j [0]);
         get_root_mod_P (c, j [0], P, verbose);
         mpz_powm_ui (j [0], j [0], 3, P);
         *no = 1;
         break;
      case CM_INVARIANT_GAMMA3:
         j = (mpz_t*) malloc (sizeof (mpz_t));
         mpz_init (j [0]);

         mpz_init (root);
         mpz_init_set_si (d_mpz, d);
         mpz_init (tmp);
         mpz_init (tmp2);

         get_root_mod_P (c, root, P, verbose);
         mpz_powm_ui (tmp, root, 2, P);
         mpz_invert (tmp2, d_mpz, P);
         mpz_mul (root, tmp, tmp2);
         mpz_add_ui (root, root, 1728);
         mpz_mod (j [0], root, P);

         *no = 1;

         mpz_clear (root);
         mpz_clear (d_mpz);
         mpz_clear (tmp);
         mpz_clear (tmp2);
         break;
      case CM_INVARIANT_WEBER:
         j = weber_cm_get_j_mod_P (c, P, no, verbose);
         break;
      case CM_INVARIANT_DOUBLEETA:
#if 0
      case CM_INVARIANT_RAMIFIED:
#endif
         j = cm_get_j_mod_P_from_modular (c, P, no, modpoldir, verbose);
         break;
      case CM_INVARIANT_SIMPLEETA:
         j = simpleeta_cm_get_j_mod_P (c, P, no, verbose);
         break;
         default: /* should not occur */
            printf ("class_cm_get_j_mod_P called for unknown class ");
            printf ("invariant\n");
            exit (1);
   }
   cm_timer_stop (clock);
   if (verbose) {
      printf ("j:");
      for (i = 0; i < *no; i++) {
         printf (" ");
         mpz_out_str (stdout, 10, j [i]);
      }
      printf ("\n");
      printf ("--- Time for j: %.1f\n", cm_timer_get (clock));
   }

   cm_class_clear (&c);

   return j;
}

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
      cm_ntl_find_root_monic (root, c.minpoly, c.minpoly_deg, P, verbose);
   else {
      mpz_init (omega);
      mpz_init (tmp);

      minpoly_P = (mpz_t*) malloc (c.minpoly_deg * sizeof (mpz_t));
      for (i = 0; i < c.minpoly_deg; i++)
         mpz_init (minpoly_P [i]);

      /* compute the second element of the integral basis modulo P */
      cm_timer_start (clock);
      fund = cm_classgroup_fundamental_discriminant (c.d);
      cm_nt_mpz_tonelli (omega, fund, P, P);
      cm_timer_stop (clock);
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

      cm_ntl_find_root_monic (root, minpoly_P, c.minpoly_deg, P, verbose);

      mpz_clear (omega);
      mpz_clear (tmp);
      for (i = 0; i < c.minpoly_deg; i++)
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

static mpz_t* cm_get_j_mod_P_from_modular (cm_class_t c, mpz_t P, int *no,
   const char* modpoldir, bool verbose)
   /* computes the possible j-values as roots of the modular polynomial      */
   /* specified by f                                                         */

{
   mpz_t  *j;
   mpz_t  omega;
      /* root of the class polynomial */
   mpz_t* poly_j;
      /* a polynomial one of whose roots is j mod P */
   int    n, i;

   mpz_init (omega);
   get_root_mod_P (c, omega, P, verbose);
   poly_j = cm_modpol_read_specialised_mod (&n, (c.p / 100) * (c.p % 100),
      CM_MODPOL_DOUBLEETA, P, omega, modpoldir);
   j = cm_ntl_find_roots (poly_j, n, P, no);
   mpz_clear (omega);
   for (i = 0; i <= n; i++)
      mpz_clear (poly_j [i]);
   free (poly_j);

   return j;
}

/*****************************************************************************/

static mpz_t* weber_cm_get_j_mod_P (cm_class_t c, mpz_t P, int *no,
   bool verbose)

{
   mpz_t* j = (mpz_t*) malloc (sizeof (mpz_t));
   mpz_t   root, f24, tmp;
   /* The following variables are only needed if the discriminant is 5 mod 8.*/
   mpz_t factor [4];
   int   i;
   mpfpx_t root_poly, f24_poly, tmp_poly, tmp2_poly;
   const unsigned long int x2 [] = {0, 0, 1};
   const unsigned long int x6 [] = {0, 0, 0, 0, 0, 0, 1};
   mpfp_t one;

   mpz_init (j [0]);
   mpz_init (root);
   mpz_init (f24);
   mpz_init (tmp);

   if (c.p != 15) {
      get_root_mod_P (c, root, P, verbose);

      if (c.d % 3 == 0)
         mpz_powm_ui (f24, root, 2ul, P);
      else
         mpz_powm_ui (f24, root, 6ul, P);

      if (c.p != 11) {
         if (c.p == 1) {
            mpz_mul_2exp (tmp, f24, 3ul);
            mpz_mod (tmp, tmp, P);
            mpz_powm_ui (f24, tmp, 4ul, P);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p == 2 || c.p == 6) {
            mpz_mul_2exp (tmp, f24, 3ul);
            mpz_mod (tmp, tmp, P);
            mpz_powm_ui (f24, tmp, 2ul, P);
            mpz_add_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p == 3) {
            mpz_mul_2exp (tmp, f24, 6ul);
            mpz_set (f24, tmp);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p == 4) {
            mpz_mul_2exp (tmp, f24, 9ul);
            mpz_set (f24, tmp);
            mpz_add_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else if (c.p == 5) {
            mpz_powm_ui (tmp, f24, 4ul, P);
            mpz_set (f24, tmp);
            mpz_sub_ui (j [0], f24, 16ul);
            mpz_powm_ui (j [0], j [0], 3ul, P);
            mpz_invert (tmp, f24, P);
            mpz_mul (j [0], j [0], tmp);
            mpz_mod (j [0], j [0], P);
         }
         else {
            mpz_mul_2exp (tmp, f24, 3ul);
            mpz_mod (tmp, tmp, P);
            mpz_powm_ui (f24, tmp, 2ul, P);
            mpz_sub_ui (j [0], f24, 16ul);
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
      for (i = 0; i <= 3; i++)
         mpz_init (factor [i]);
      mpfpx_type_init (P, MPFPX_KARATSUBA);
      mpfpx_init (root_poly);
      mpfpx_init (f24_poly);
      mpfpx_init (tmp_poly);
      mpfpx_init (tmp2_poly);
      mpfp_init_set_ui (one, 1ul);

      cm_ntl_find_factor (factor, c.minpoly, c.minpoly_deg, 3, P, verbose);
      mpfpx_set_z_array (root_poly, factor, 4);
      printf ("Root: "); mpfpx_out (root_poly);

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

   mpz_clear (root);
   mpz_clear (f24);
   mpz_clear (tmp);

   return j;
}

/*****************************************************************************/

static mpz_t* simpleeta_cm_get_j_mod_P (cm_class_t c, mpz_t P, int *no,
   bool verbose)

{
   mpz_t* j = (mpz_t*) malloc (sizeof (mpz_t));
   mpz_t   root, f3, tmp;
   int l = c.p / 10000;

   mpz_init (j [0]);
   mpz_init (root);
   mpz_init (f3);
   mpz_init (tmp);

   get_root_mod_P (c, root, P, verbose);
   /* raise to the power s/e */
   mpz_powm_ui (root, root,
      (unsigned long int) ((c.p % 100) / ((c.p / 100) % 100)), P);

   if (l == 3) {
      mpz_add_ui (f3, root, 3ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_invert (tmp, root, P);
      mpz_mul_ui (tmp, tmp, 27ul);
      mpz_add_ui (tmp, tmp, 1ul);
      mpz_mul (j [0], f3, tmp);
      mpz_mod (j [0], j [0], P);
   }
   else if (l == 5) {
      mpz_add_ui (tmp, root, 10ul);
      mpz_mul (f3, root, tmp);
      mpz_mod (f3, f3, P);
      mpz_add_ui (f3, f3, 5ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_invert (tmp, root, P);
      mpz_mul (j [0], f3, tmp);
      mpz_mod (j [0], j [0], P);
   }
   else if (l == 7) {
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
   else if (l == 13) {
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
   else if (l == 4) {
      mpz_add_ui (f3, root, 16ul);
      mpz_mul (f3, f3, root);
      mpz_mod (f3, f3, P);
      mpz_invert (tmp, f3, P);
      mpz_add_ui (f3, f3, 16ul);
      mpz_powm_ui (f3, f3, 3ul, P);
      mpz_mul (j [0], tmp, f3);
      mpz_mod (j [0], j [0], P);
   }
   else if (l == 9) {
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
   else if (l == 25) {
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
      mpz_mod (j [0], P, P);
   }

   *no = 1;

   mpz_clear (root);
   mpz_clear (f3);
   mpz_clear (tmp);

   return j;
}

/*****************************************************************************/
/*****************************************************************************/
