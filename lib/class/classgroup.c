/*

classgroup.c - computations with class groups and quadratic forms

Copyright (C) 2009, 2010, 2018 Andreas Enge

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

typedef struct avl_t {
   struct avl_t *l, *r; /* left and right subtrees */
   signed char     b; /* balance factor; -1, 0 or 1 */
   cm_form_t       c; /* content */
} avl_t;


static int_cl_t classgroup_gcdext (int_cl_t *u, int_cl_t *v, int_cl_t a,
   int_cl_t b);
static int_cl_t classgroup_tonelli (int_cl_t a, int_cl_t p);
static int_cl_t classgroup_fundamental_discriminant_conductor (int_cl_t d,
   uint_cl_t *cond_primes, unsigned int *cond_exp);
static int avl_cmp (cm_form_t P, cm_form_t Q);
static bool avl_insert_rec (avl_t **t, cm_form_t c, bool *ok);
static bool avl_insert (avl_t **t, cm_form_t c);
static int avl_flatten_rec (cm_form_t *list, avl_t *t, int pos);
static void avl_flatten (cm_form_t *list, avl_t *t);
static void avl_delete (avl_t *t);

/*****************************************************************************/
/*                                                                           */
/* functions transforming between int_cl_t and mpz_t                         */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_mpz_set_icl (mpz_t rop, int_cl_t op)
   /* relies on the fact that int_cl_t has 64 bits */

{
   int_cl_t mask = (((int_cl_t) 1) << 32) - 1;
   int_cl_t abso = (op > 0 ? op : -op);

   /* copy 32 high bits */
   mpz_set_ui (rop, (unsigned long int) (abso >> 32));
   /* add 32 low bits */
   mpz_mul_2exp (rop, rop, 32);
   mpz_add_ui (rop, rop, (unsigned long int) (abso & mask));

   if (op < 0)
      mpz_neg (rop, rop);
}

/*****************************************************************************/

int_cl_t cm_classgroup_mpz_get_icl (mpz_t op)

{
   int_cl_t rop;
   int sign = (mpz_cmp_ui (op, 0) < 0 ? -1 : 1);
   mpz_t tmp;

   mpz_init (tmp);

   /* switch to absolute value */
   if (sign < 0)
      mpz_neg (op, op);

   /* extract 32 high bits */
   mpz_tdiv_q_2exp (tmp, op, 32);
   rop = ((int_cl_t) (mpz_get_ui (tmp))) << 32;
   /* extract 32 low bits */
   mpz_tdiv_r_2exp (tmp, op, 32);
   rop += mpz_get_ui (tmp);

   if (sign < 0) {
      rop = -rop;
      mpz_neg (op, op);
   }

   mpz_clear (tmp);

   return rop;
}


/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_classgroup_init (cm_classgroup_t *cl, int_cl_t disc, bool verbose)

{
   int h;
   cm_form_t *Cl;
      /* class number and class group as far as they are already computed */
   avl_t *t;
      /* avl tree of forms in Cl */
   int_cl_t p;
   cm_form_t P, Ppow;
      /* prime, form above it and powers of the form */
   int relo;
      /* relative order of the prime form modulo Cl */
   int i, j;

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

   cl->h = cm_classgroup_h (cl->d);
   cl->form = (cm_form_t *) malloc (cl->h * sizeof (cm_form_t));
   cl->conj = (int *) malloc (cl->h * sizeof (int));
   if (verbose)
      printf ("Class number: h = %d\n", cl->h);

   /* computation of the class group */
   Cl = (cm_form_t *) malloc (cl->h * sizeof (cm_form_t));
   P.a = 1;
   P.b = (disc % 2 == 0 ? 0 : 1);
   t = NULL;
   avl_insert (&t, P);
   h = 1;
   Cl [0] = P;

   p = (int_cl_t) 1;
   while (h < cl->h) {
      /* select next prime form */
      do {
         p = (int_cl_t) cm_nt_next_prime ((unsigned long int) p);
      } while (cm_classgroup_kronecker (disc, p) == -1 ||
           /* avoid primes dividing the conductor */
           (cm_classgroup_mod (disc, p*p) == 0
            && (   cm_classgroup_mod (disc / (p*p), (uint_cl_t) 4) == 1
                || cm_classgroup_mod (disc / (p*p), (uint_cl_t) 4) == 0)));
      P = cm_classgroup_prime_form (p, disc);

      /* determine its relative order through inserting its powers into   */
      /* the tree                                                         */
      Ppow = P;
      relo = 1;
      while (avl_insert (&t, Ppow)) {
         cm_classgroup_compose (&Ppow, Ppow, P, disc);
         relo++;
      }
      if (verbose && relo > 1)
         printf ("   [%"PRIicl", %"PRIicl"]: %i\n", P.a, P.b, relo);

      /* Multiply all other forms by the powers of P and insert them into
         the tree. */
      for (i = 0; i < h; i++) {
         Ppow = Cl [i];
         for (j = 1; j < relo; j++) {
            cm_classgroup_compose (&Ppow, Ppow, P, disc);
            avl_insert (&t, Ppow);
         }
      }

      h *= relo;
      /* copy tree content into array */
      avl_flatten (Cl, t);
   }

   for (i = 0; i < cl->h; i++)
      cl->form [i] = Cl [i];

   /* Pair up inverse forms. */
   for (i = 0; i < cl->h; i++) {
      for (j = i;
           j < cl->h && cl->form [j].b != -cl->form [i].b;
           j++);
      if (j == cl->h || cl->form [i].b == 0)
         cl->conj [i] = i;
      else {
         cl->conj [i] = j;
         cl->conj [j] = i;
      }
   }

   avl_delete (t);
   free (Cl);
}

/*****************************************************************************/

void cm_classgroup_clear (cm_classgroup_t *cl)

{
   free (cl->form);
   free (cl->conj);
}

/*****************************************************************************/
/*                                                                           */
/* elementary number theory with int_cl_t                                    */
/*                                                                           */
/*****************************************************************************/

uint_cl_t cm_classgroup_mod (int_cl_t a, uint_cl_t p)
      /* returns a representative of a % p in [0, p-1[ */

{
   int_cl_t res = a % (int_cl_t) p;
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
   /* in Cohen93 (by the binary algorithm)                                   */
{
   const int tab [8] = {0, 1, 0, -1, 0, -1, 0, 1};
   int k;
   int_cl_t r;

   /* step 1 */
   if (b == 0) {
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
   if (b % 2 == 0) {
      b >>= 1;
      k = tab [cm_classgroup_mod (a, (uint_cl_t) 8)];
   }
   else
      k = 1;

   if (b < 0) {
      b = -b;
      if (a < 0)
         k = -k;
   }

   /* step 3 */
   a = cm_classgroup_mod (a, b);

   /* steps 4 to 6 */
   while (a != 0) {
      while (a % 4 == 0)
         a >>= 2;
      if (a % 2 == 0) {
         a >>= 1;
         k *= tab [cm_classgroup_mod (b, (uint_cl_t) 8)];
      }
      if (b > a) {
         r = b - a;
         if (   cm_classgroup_mod (a, (uint_cl_t) 4) == 3
             && cm_classgroup_mod (b, (uint_cl_t) 4) == 3)
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

static int_cl_t classgroup_tonelli (int_cl_t a, int_cl_t p)
{
   mpz_t az, pz, rz;
   int_cl_t r;

   mpz_init (az);
   mpz_init (pz);
   mpz_init (rz);

   cm_classgroup_mpz_set_icl (az, a);
   cm_classgroup_mpz_set_icl (pz, p);
   cm_nt_mpz_tonelli_z (rz, az, pz);
   r = cm_classgroup_mpz_get_icl (rz);

   mpz_clear (az);
   mpz_clear (pz);
   mpz_clear (rz);

   return r;
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
   int pow4;
   uint_cl_t factors [17];
   unsigned int exponents [17];
   int_cl_t fundamental_d = -1;

   /* handle 2 in the conductor separately */
   local_d = d;
   pow4 = 0;
   while (local_d % 4 == 0) {
      local_d /= 4;
      if ((local_d - 1) % 4 == 0 || local_d % 4 == 0)
         pow4++;
      else
         fundamental_d *= 4;
   }
   if (local_d % 2 == 0) {
      local_d /= 2;
      fundamental_d *= 2;
   }

   if (pow4 != 0) {
      cond_primes [0] = 2;
      cond_exp [0] = pow4;
      j = 1;
   }
   else
      j = 0;

   cm_classgroup_factor (local_d, factors, exponents);

   for (i = 0; factors [i] != 0; i++) {
      if (exponents [i] >= 2) {
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

int cm_classgroup_h (int_cl_t d)
   /* computes the class number of the imaginary quadratic order with        */
   /* discriminant d using Louboutin's formula [Lou02] for the maximal part  */
   /* and the class number formula for non-maximal orders.                   */

{
   int_cl_t fund;
   uint_cl_t cond_primes [17];
   unsigned int cond_exp [17];
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

   return h;
}

/*****************************************************************************/
/*                                                                           */
/* functions implementing avl trees on classgroup elements                   */
/*                                                                           */
/*****************************************************************************/

static int avl_cmp (cm_form_t P, cm_form_t Q)
   /* Return -1, 0 or 1, depending on whether P is smaller than, equal to
      or larger than Q. Uses the lexicographical order on (a, b). */

{
   if (P.a < Q.a)
      return 1;
   else if (P.a > Q.a)
      return -1;
   else {
      if (P.b < Q.b)
         return -1;
      else if (P.b > Q.b)
         return 1;
      else
         return 0;
   }
}

/*****************************************************************************/

static bool avl_insert_rec (avl_t **t, cm_form_t c, bool *ok)
   /* Recursive function doing the work; the return value indicates whether  */
   /* the branch has become longer, so that rebalancing is needed.           */

{
   bool rebalance;
   int cmp;
   avl_t *y, *x, *w;
      /* notations as in the GNU libavl book by Ben Pfaff */

   y = *t;
   if (y == NULL) {
      *t = (avl_t *) malloc (sizeof (avl_t));
      (*t)->l = NULL;
      (*t)->r = NULL;
      (*t)->b = 0;
      (*t)->c = c;
      *ok = true;
      return true;
   }

   cmp = avl_cmp (c, y->c);
   if (cmp == 0) {
      *ok = false;
      return false;
   }
   else if (cmp < 0) {
      /* insert into the left subtree */
      rebalance = avl_insert_rec (&(y->l), c, ok);
      if (rebalance) {
         /* left subtree has become longer */
         switch (y->b) {
            case 1: /* right branch was longer, balanced now */
               y->b = 0;
               return false;
            case 0: /* now left branch is longer, rebalance needed */
               y->b = -1;
               return true;
            default: /* now left branch is longer by 2, rotations needed */
               x = y->l;
               if (x->b == -1) {
                  /* rotate right at y */
                  y->l = x->r;
                  x->r = y;
                  *t = x;
                  /* update balance factor */
                  y->b = 0;
                  x->b = 0;
               }
               else {
                  /* x->b == 1; rotate left at x, then right at y */
                  w = x->r;
                  x->r = w->l;
                  y->l = w->r;
                  w->l = x;
                  w->r = y;
                  *t = w;
                  if (w->b == 1)
                     x->b = -1;
                  else
                     x->b = 0;
                  if (w->b == -1)
                     y->b = 1;
                  else
                     y->b = 0;
                  w->b = 0;
               }
               /* everything is balanced now */
               return false;
         }
      }
      else /* left subtree has same length as before */
         return false;
   }
   else {
      /* c > y->c, insert into the right subtree and balance by symmetry */
      rebalance = avl_insert_rec (&(y->r), c, ok);
      if (rebalance) {
         switch (y->b) {
            case -1:
               y->b = 0;
               return false;
            case 0:
               y->b = 1;
               return true;
            default:
               x = y->r;
               if (x->b == 1) {
                  y->r = x->l;
                  x->l = y;
                  *t = x;
                  y->b = 0;
                  x->b = 0;
               }
               else {
                  w = x->l;
                  x->l = w->r;
                  y->r = w->l;
                  w->r = x;
                  w->l = y;
                  *t = w;
                  if (w->b == -1)
                     x->b = 1;
                  else
                     x->b = 0;
                  if (w->b == 1)
                     y->b = -1;
                  else
                     y->b = 0;
                  w->b = 0;
               }
               return false;
         }
      }
      else
         return false;
   }
}

/*****************************************************************************/

static bool avl_insert (avl_t **t, cm_form_t c)
   /* Tries to insert a new element c into the tree t; returns FALSE if the  */
   /* element is already contained in the tree, TRUE otherwise (indicating   */
   /* that an insertion has indeed taken place).                             */

{
   bool ok;
   avl_insert_rec (t, c, &ok);
   return ok;
}

/*****************************************************************************/

static int avl_flatten_rec (cm_form_t *list, avl_t *t, int pos)
   /* Returns the content of t as an ordered array from position pos in      */
   /* list, which needs to provide the necessary space. The return value is  */
   /* the next free index.                                                   */

{
   if (t == NULL)
      return pos;
   else {
      pos = avl_flatten_rec (list, t->l, pos);
      list [pos] = t->c;
      pos++;
      return avl_flatten_rec (list, t->r, pos);
   }
}

/*****************************************************************************/

static void avl_flatten (cm_form_t *list, avl_t *t)
   /* Returns the content of t as an ordered array in list, which needs to   */
   /* provide the necessary space.                                           */

{
   avl_flatten_rec (list, t, 0);
}

/*****************************************************************************/

static void avl_delete (avl_t *t)

{
   if (t != NULL) {
      avl_delete (t->l);
      avl_delete (t->r);
      free (t);
   }
}

/*****************************************************************************/
/*                                                                           */
/* functions for computing in the class group                                */
/*                                                                           */
/*****************************************************************************/

int_cl_t cm_classgroup_compute_c (int_cl_t a, int_cl_t b, int_cl_t d)
   /* computes c = (b^2 - d) / (4*a), potentially switching to multi-        */
   /* precision to avoid intermediate overflow                               */

{
   int_cl_t c;

   if ((b > 0 ? b : -b)
         < ((int_cl_t) 1) << (4 * sizeof (int_cl_t) - 2))
      c = (b * b - d) / a / 4;
   else {
      /* should happen only rarely */
      mpz_t az, dz, cz;
      mpz_init (az);
      mpz_init (cz);
      mpz_init (dz);
      cm_classgroup_mpz_set_icl (dz, d);
      cm_classgroup_mpz_set_icl (az, a);
      cm_classgroup_mpz_set_icl (cz, b);
      mpz_mul (cz, cz, cz);
      mpz_sub (cz, cz, dz);
      mpz_div (cz, cz, az);
      mpz_div_2exp (cz, cz, 2);
      c = cm_classgroup_mpz_get_icl (cz);
      mpz_clear (az);
      mpz_clear (cz);
      mpz_clear (dz);
   }

   return c;
}

/*****************************************************************************/

void cm_classgroup_reduce (cm_form_t *Q, int_cl_t d)
   /* reduces the quadratic form Q without checking if it belongs indeed to */
   /* the discriminant d.                                                   */

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

      c = cm_classgroup_compute_c (Q->a, Q->b, d);

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
   /* Q2.                                                                     */

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

cm_form_t cm_classgroup_prime_form (int_cl_t p, int_cl_t d)
   /* Assumes that p is a ramified or split prime and returns the reduction  */
   /* of the prime form above p with non-negative b.                         */

{
   cm_form_t Q;

   Q.a = p;
   if (p == 2)
      if (d % 8 == 0)
         Q.b = 0;
      else if ((d - 4) % 8 == 0)
         Q.b = 2;
      else
         Q.b = 1;
   else {
      Q.b = classgroup_tonelli (d, p);
      /* fix parity of b */
      if ((d + Q.b) % 2 != 0)
         Q.b += p;
   }

   cm_classgroup_reduce (&Q, d);

   return Q;
}

/*****************************************************************************/
/*****************************************************************************/
