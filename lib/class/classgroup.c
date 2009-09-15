#include "cm_class-impl.h"

#define mpz_add_si(c, a, b) \
   (b >= 0) ? mpz_add_ui (c, a, b) \
            : mpz_sub_ui (c, a, (unsigned long int) (-b))
#define mpz_sub_si(c, a, b) \
   (b >= 0) ? mpz_sub_ui (c, a, b) \
            : mpz_add_ui (c, a, (unsigned long int) (-b))


static int_cl_t classgroup_gcdext (int_cl_t *u, int_cl_t *v, int_cl_t a,
   int_cl_t b);
static int_cl_t classgroup_fundamental_discriminant_conductor (int_cl_t d,
   uint_cl_t *cond_primes, unsigned int *cond_exp);
static void classgroup_write (cm_classgroup_t cl);
static bool classgroup_read (cm_classgroup_t cl);

/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_init (cm_classgroup_t *cl, int_cl_t disc, bool checkpoints,
   bool verbose)
   /* If checkpoints is true, then the function tries to read the         */
   /* class group from a file and writes the result to a file             */

{
   int k;
   int_cl_t a, b, c;

   if (disc >= 0) {
      printf ("\n*** The discriminant must be negative.\n");
      exit (1);
   }
   else if (disc % 4 != 0 && (disc - 1) % 4 != 0) {
      printf ("\n*** %"PRIicl" is not a quadratic discriminant.\n", disc);
      exit (1);
   }
   else
      cl->d = disc;

   cl->h = cm_classgroup_h (&(cl->h1), &(cl->h2), cl->d);
   cl->h12 = cl->h1 + cl->h2;
   cl->form = (cm_form_t *) malloc (cl->h12 * sizeof (cm_form_t));

   if (!checkpoints || !classgroup_read (*cl)) {
      k = 0;
      for (a = 1; a*a <= (-cl->d) / 3; a++) {
         if (cl->d % 2 == 0)
            b = 0;
         else
            b = 1;
         for (; b <= a; b += 2) {
            if ((b*b - cl->d) % (4*a) == 0) {
               c = ((b*b - cl->d) / (4*a));
               if (c >= a) {
                  if (cm_classgroup_gcd (cm_classgroup_gcd (a, b), c) == 1) {
                     /* we have a primitive reduced form */
                     cl->form [k].a = a;
                     cl->form [k].b = b;
                     if (b == 0 || b == a || a == c)
                        cl->form [k].emb = real;
                     else
                        cl->form [k].emb = complex;
                     k++;
                  }
               }
            }
         }
      }

   if (checkpoints)
      classgroup_write (*cl);
   }

   if (verbose)
      printf ("Class numbers: h = %d, h1 = %d, h2 = %d\n",
         cl->h, cl->h1, cl->h2);
}

/*****************************************************************************/

void cm_classgroup_clear (cm_classgroup_t *cl)

{
   free (cl->form);
}

/*****************************************************************************/
/*                                                                           */
/* file handling functions                                                   */
/*                                                                           */
/*****************************************************************************/

static void classgroup_write (cm_classgroup_t cl)
   /* writes the classgroup to the file                                      */
   /* CLASS_TMPDIR + "/tmp_" + d + "_classgroup.dat"                         */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_classgroup.dat", CM_CLASS_TMPDIR, -cl.d);

   if (!cm_file_open_write (&f, filename))
      exit (1);

   fprintf (f, "%"PRIicl"\n", -cl.d);
   fprintf (f, "%i\n", cl.h);
   fprintf (f, "%i\n", cl.h1);
   fprintf (f, "%i\n", cl.h2);

   for (i = 0; i < cl.h12; i++)
      fprintf (f, "%"PRIicl" %"PRIicl"\n",
         cl.form [i].a, cl.form [i].b);

   cm_file_close (f);
}

/*****************************************************************************/

bool classgroup_read (cm_classgroup_t cl)
   /* reads the classgroup from the file                                     */
   /* CLASS_TMPDIR + "/tmp_" + d + "_classgroup.dat"                         */
   /* If an error occurs, the return value is false.                         */

{
   char filename [255];
   FILE* f;
   int i;
   int_cl_t disc;

   sprintf (filename, "%s/tmp_%"PRIicl"_classgroup.dat", CM_CLASS_TMPDIR, -cl.d);

   if (!cm_file_open_read (&f, filename))
      return false;

   if (!fscanf (f, "%"SCNicl"\n", &disc))
      return false;
   if (-disc != cl.d) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** discriminant %"PRIicl" instead of %"PRIicl"\n", -disc, cl.d);
      exit (1);
   }
   if (!fscanf (f, "%i", &i))
      return false;
   if (i != cl.h) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** h equals %i instead of %i\n", i, cl.h);
      exit (1);
   }
   if (!fscanf (f, "%i", &i))
      return false;
   if (i != cl.h1) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** h1 equals %i instead of %i\n", i, cl.h1);
      exit (1);
   }
   if (!fscanf (f, "%i", &i))
      return false;
   if (i != cl.h2) {
      printf ("*** Inconsistency between file '%s' ", filename);
      printf ("and internal data:\n");
      printf ("*** h2 equals %i instead of %i\n", i, cl.h2);
      exit (1);
   }

   for (i = 0; i < cl.h12; i++) {
      if (!fscanf (f, "%"SCNicl" %"SCNicl,
           &(cl.form [i].a), &(cl.form [i].b)))
         return false;
      if (cl.form [i].b == 0 || cl.form [i].b == cl.form [i].a
          || cl.form [i].a == ((cl.form [i].b*cl.form [i].b - cl.d) / (4*cl.form [i].a)))
         cl.form [i].emb = real;
      else
         cl.form [i].emb = complex;
   }

   cm_file_close (f);

   return true;
}

/*****************************************************************************/
/*                                                                           */
/* elementary number theory with int_cl_t                                    */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_mpz_set_icl (mpz_t rop, int_cl_t op)

{
   char op_str [22];
      /* 64 bits = 20 decimal digits, plus sign, plus '\0' */

   sprintf (op_str, "%"PRIicl, op);
   mpz_set_str (rop, op_str, 10);
}

/*****************************************************************************/

uint_cl_t cm_classgroup_mod (int_cl_t a, uint_cl_t p)
      /* returns a representative of a % p in [0, p-1[ */

{
   int_cl_t res = a % p;
   return (res < 0 ? res + (int_cl_t) p : res);
}

/*****************************************************************************/

int_cl_t cm_classgroup_gcd (int_cl_t a, int_cl_t b)
   /* returns the positive gcd of a and b */

{
   int_cl_t r;

   if (a == 0)
      return b;
   else if (b == 0)
      return a;
   else
   {
      if (a < 0)
         a = -a;
      if (b < 0)
         b = -b;
      r = a % b;
      while (r > 0)
      {
         a = b;
         b = r;
         r = a % b;
      }
      return b;
   }
}

/*****************************************************************************/

static int_cl_t classgroup_gcdext (int_cl_t *u, int_cl_t *v, int_cl_t a,
   int_cl_t b)
   /* returns the positive gcd d of a and b; if u and v is not NULL,         */
   /* modifies them such that u a + v b = d; it is also possible to have     */
   /* only one of u and v non NULL.                                          */

{
   int_cl_t r0, r1, r2, u0, u1, u2, v0, v1, v2, q;
   int sgn_a, sgn_b;

   if (a < 0) {
      sgn_a = -1;
      r0 = -a;
   }
   else {
      sgn_a = 1;
      r0 = a;
   }
   if (b < 0) {
      sgn_b = -1;
      r1 = -b;
   }
   else {
      sgn_b = 1;
      r1 = b;
   }
   /* computing u and v might be faster than permanent tests for NULL */
   u0 = 1;
   u1 = 0;
   v0 = 0;
   v1 = 1;

   while (r1 != 0) {
      q = r0 / r1;
      r2 = r0 % r1;
      r0 = r1;
      r1 = r2;
      u2 = u0 - q * u1;
      u0 = u1;
      u1 = u2;
      v2 = v0 - q * v1;
      v0 = v1;
      v1 = v2;
   }

   if (u != NULL)
      *u = sgn_a * u0;
   if (v != NULL)
      *v = sgn_b * v0;

   return r0;
}

/*****************************************************************************/

int cm_classgroup_kronecker (int_cl_t a, int_cl_t b)
   /* computes the Kronecker symbol (a / b) following Algorithm 1.4.12       */
   /* in Cohen (by the binary algorithm)                                     */
{
   const int tab [8] = {0, 1, 0, -1, 0, -1, 0, 1};
   int k;
   int_cl_t r;

   /* step 1 */
   if (b == 0)
   {
      if (a == 1 || a == -1)
         return 1;
      else
         return 0;
   }

   /* step 2*/
   if (a % 2 == 0 && b % 2 == 0)
      return 0;

   while (b % 4 == 0)
      b >>= 2;
   if (b % 2 == 0)
   {
      b >>= 1;
      k = tab [a & 7];
   }
   else
      k = 1;

   if (b < 0)
   {
      b = -b;
      if (a < 0)
         k = -k;
   }

   /* step 3 and added test; here, b is already odd */
   if (a < 0)
   {
      a = -a;
      if (b & 2)
         k = -k;
   }
   a %= b;

   /* steps 4 to 6 */
   while (a != 0)
   {
      while (a % 4 == 0)
         a >>= 2;
      if (a % 2 == 0)
      {
         a >>= 1;
         k *= tab [b & 7];
      }
      if (b > a)
      {
         r = b - a;
         if (a & b & 2)
            k = -k;
         b = a;
         a = r;
      }
      else
         a -= b;
   }

   if (b > 1)
      return 0;
   else
      return k;
}

/*****************************************************************************/
/*                                                                           */
/* functions computing fundamental discriminants, conductors and class       */
/* numbers                                                                   */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_factor (int_cl_t d, uint_cl_t *factors,
                   unsigned int *exponents)
   /* factors the absolute value of d by trial division. The prime factors   */
   /* are stored in "factors", their multiplicities in "exponents", which    */
   /* must provide sufficiently much space. The list of prime factors is     */
   /* terminated by 0, so that 12 entries suffice for a number of 32 bits,   */
   /* and 17 entries for a number of 64 bits.                                */

{
   uint_cl_t no, trial, trial2;
   int j;

   if (d < 0)
      no = -d;
   else
      no = d;

   j = 0;
   trial = 0;
   trial2 = 0;
   while (trial2 <= no) {
      if (trial == 0) {
         trial = 2;
         trial2 = 4;
      }
      else if (trial == 2) {
         trial = 3;
         trial2 = 9;
      }
      else {
         trial += 2;
         trial2 += 4 * (trial - 1);
      }
      if (no % trial == 0) {
         factors [j] = trial;
         no /= trial;
         exponents [j] = 1;
         while (no % trial == 0) {
            no /= trial;
            exponents [j]++;
         }
         j++;
      }
   }
   if (no != 1) {
     factors [j] = no;
     exponents [j] = 1;
     j++;
   }
   factors [j] = 0;
}

/*****************************************************************************/

static int_cl_t classgroup_fundamental_discriminant_conductor (int_cl_t d,
   uint_cl_t *cond_primes, unsigned int *cond_exp)
   /* returns the fundamental discriminant corresponding to d, and the       */
   /* prime factorisation of the conductor via cond_primes and cond_exp.     */

{
   int i, j;
   int_cl_t local_d;
   int exp2;
   uint_cl_t factors [17];
   unsigned int exponents [17];
   int_cl_t fundamental_d = -1;

   /* handle 2 in the conductor separately */
   local_d = d;
   exp2 = 0;
   while (local_d % 4 == 0)
   {
      local_d /= 4;
      if ((local_d - 1) % 4 == 0 || local_d % 4 == 0)
         exp2++;
      else
         fundamental_d *= 4;
   }
   if (local_d % 2 == 0)
   {
      local_d /= 2;
      fundamental_d *= 2;
   }

   if (exp2 != 0)
   {
      cond_primes [0] = 2;
      cond_exp [0] = exp2;
      j = 1;
   }
   else
      j = 0;

   cm_classgroup_factor (local_d, factors, exponents);

   for (i = 0; factors [i] != 0; i++)
   {
      if (exponents [i] >= 2)
      {
         cond_primes [j] = factors [i];
         cond_exp [j] = exponents [i] / 2;
         j++;
      }
      if (exponents [i] & 1)
            fundamental_d *= factors [i];
   }
   cond_primes [j] = 0;

   return fundamental_d;
}

/*****************************************************************************/

int_cl_t cm_classgroup_fundamental_discriminant (int_cl_t d)

{
   uint_cl_t cond_primes [17];
   unsigned int cond_exp [17];

   return classgroup_fundamental_discriminant_conductor (d, cond_primes,
             cond_exp);
}

/*****************************************************************************/

int cm_classgroup_h (int *h1, int *h2, int_cl_t d)
   /* computes the class number of the imaginary quadratic order with        */
   /* discriminant d using Louboutin's formula [Lou02] for the maximal part  */
   /* and the class number formula for non-maximal orders.                   */
   /* If h1 is not the NULL pointer, the numbers of ambiguous reduced forms  */
   /* and of pairs of non-ambiguous ones are returned via h1 and h2.         */

{
   int_cl_t fund;
   uint_cl_t cond_primes [17];
   unsigned int cond_exp [17];
   uint_cl_t factors [17];
   double   pi, a2, a, en, enp1, fn, deltaf, sum1, sum2;
      /* en stands for e_n = exp (- n^2 a2), enp1 for e_{n+1},               */
      /* deltaf for exp (-2 a2) and  fn for exp (- (2n+3) a2),               */
      /* so that e_{n+2} = e_{n+1} * fn.                                     */
      /* sum1 contains the sum of chi (n) * e_n / n;                         */
      /* sum2 contains the sum of (e_n + e_{n+1}) S_n.                       */
   int m, n, S, h;
   int chi, i;

   fund = classgroup_fundamental_discriminant_conductor (d, cond_primes,
      cond_exp);

   if (fund == -3 || fund == -4)
      h = 1;
   else
   {
      pi = 2 * asin (1.0);
      a2 = pi / (-fund);
      a = sqrt (a2);
      m = (int) ceil (sqrt (-fund / (2 * pi) * log (-fund / pi) + 6));

      /* initialisation for n = 1 */
      chi = 1;
      S = 1;
      en = exp (-a2);
      deltaf = en * en;
      enp1 = deltaf * deltaf;
      fn = enp1 * en;
      sum1 = en;
      sum2 = en + enp1;

      for (n = 2; n <= m; n++)
      {
         chi = cm_classgroup_kronecker (fund, (int_cl_t) n);
         S += chi;
         en = enp1;
         enp1 *= fn;
         fn *= deltaf;
         sum1 += chi * en / (double) n;
         sum2 += (en + enp1) * S;
      }

      h = (int) ((sum1 / a + sum2 * a) / sqrt (pi) + 0.5);
   }

   /* correct for the conductor */
   for (i = 0; cond_primes [i] != 0; i++)
   {
      h *=   cond_primes [i]
           - cm_classgroup_kronecker (fund, (int_cl_t) cond_primes [i]);
      while (cond_exp [i] > 1)
      {
         h *= cond_primes [i];
         cond_exp [i]--;
      }
   }

   /* correct for the units */
   if (cond_primes [0] != 0) {
      if (fund == -3)
         h /= 3;
      else if (fund == -4)
         h /= 2;
   }

   if (h1 != NULL) {
      /* compute the number of ambiguous forms; factors the discriminant     */
      /* again, which is not very elegant, but fast enough                   */
      m = 0;
      fund = d;
      while (fund % 4 == 0) {
         fund /= 4;
         m++;
      }
      if ((fund - 1) % 4 == 0) {
         if (m <= 1)
            *h1 = 1;
         else if (m == 2)
            *h1 = 2;
         else
            *h1 = 4;
      }
      else if ((fund - 3) % 4 == 0) {
         if (m <= 2)
            *h1 = 2;
         else
            *h1 = 4;
      }
      else {
         if (m <= 1)
            *h1 = 1;
         else
            *h1 = 2;
      }

      if (fund != -1) {
         cm_classgroup_factor (fund, factors, cond_exp);
         for (i = 1; factors [i] != 0; i++)
            *h1 *= 2;
      }
      else
         *h1 /= 2;

      *h2 = (h - *h1) / 2;
   }

   return h;
}

/*****************************************************************************/
/*                                                                           */
/* functions for computing in the class group                                */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_reduce (cm_form_t *Q, int_cl_t d)
   /* reduces the quadratic form Q without checking if it belongs indeed to */
   /* the discriminant d and without computing Q.emb.                       */
   /* The function is of rather limited use, since it can only be used      */
   /* when Q.b is noticeably smaller than 2^32; otherwise there is an       */
   /* overflow in the computation of Q.b^2 - d.                             */

{
   int_cl_t c, a_minus_b, two_a, offset;
   bool reduced;

   reduced = false;
   while (!reduced){
      /* first step: obtain |b| <= a */
      a_minus_b = Q->a - Q->b;
      two_a = 2 * Q->a;
      if (a_minus_b < 0) {
         a_minus_b++;
         /* a trick to obtain the correct rounding */
         offset = a_minus_b / two_a;
         /* Since a_minus_b is negative, a negative remainder is computed, */
         /* so the quotient is effectively rounded to the nearest integer  */
         offset--;
         /* offset is (a-b) / (2a) floored */
         offset *= two_a;
         Q->b += offset;
      }
      else if (a_minus_b >= two_a) {
         offset = a_minus_b / two_a;
         offset *= two_a;
         Q->b += offset;
      }

      assert (Q->b < ((int_cl_t) 1) << (4 * sizeof (int_cl_t) - 2));
         /* prevent overflow in the computation of c = b^2 - d */
      /* compute c */
      c = (Q->b * Q->b - d) / (4 * Q->a);
      /* if not reduced, invert */
      if (Q->a < c || (Q->a == c && Q->b >= 0))
            reduced = true;
      else {
         Q->a = c;
         Q->b = -Q->b;
      }
   }
}

/*****************************************************************************/

void cm_classgroup_compose (cm_form_t *Q, cm_form_t Q1, cm_form_t Q2,
   int_cl_t d)
   /* computes the reduced form Q corresponding to the composition of Q1 and */
   /* without computing Q.emb                                                */

{
   int_cl_t s, t, v1, v, w, a2t;

   t = classgroup_gcdext (&v1, NULL, Q2.a, Q1.a);

   if (t == 1) {
      Q->a = Q1.a * Q2.a;
      Q->b = Q2.b + Q2.a * v1 * (Q1.b - Q2.b);
   }
   else {
      s = (Q1.b + Q2.b) / 2;
      t = classgroup_gcdext (&w, &v, s, t);
      v *= v1;
      a2t = Q2.a / t;
      Q->a = (Q1.a / t) * a2t;
      Q->b = ((s - Q2.b) * v - w * (Q2.b * Q2.b - d) / (4 * Q2.a))
             % (2 * Q->a); /* intermediate reduction to avoid overflow */
      Q->b = (Q2.b + 2 * Q->b * a2t) % (2 * Q->a);
   }
   cm_classgroup_reduce (Q, d);
}

/*****************************************************************************/
