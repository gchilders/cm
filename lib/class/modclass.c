#include "cm_class-impl.h"

#define mpz_sub_si(c, a, b) \
   (b >= 0) ? mpz_sub_ui (c, a, (unsigned long int) b) \
            : mpz_add_ui (c, a, (unsigned long int) (-b))

static void write_q24eta (cm_modclass_t mc, mpc_t *q24eta, char *type);
static bool read_q24eta (cm_modclass_t mc, mpc_t *q24eta, char *type);
static void compute_q24 (cm_modclass_t mc, mpc_t *eta, bool verbose);

/*****************************************************************************/
/*                                                                           */
/* creating and freeing variables                                            */
/*                                                                           */
/*****************************************************************************/

void cm_modclass_init (cm_modclass_t *mc, cm_classgroup_t cl,
   mp_prec_t prec, bool checkpoints, bool verbose)

{
   int i;
   mpc_t *q24;
   cm_timer  clock1, clock2;

   cm_modular_init (&(mc->m), prec);
   mc->cl = cl;

   cm_timer_start (clock1);
   mpfr_init2 (mc->root, prec);
   mpfr_sqrt_ui (mc->root, -mc->cl.d, GMP_RNDN);
   mpfr_init2 (mc->sqrt2_over2, prec);
   mpfr_sqrt_ui (mc->sqrt2_over2, 2ul, GMP_RNDN);
   mpfr_div_2ui (mc->sqrt2_over2, mc->sqrt2_over2, 1ul, GMP_RNDN);
   mpfr_init2 (mc->sqrt2_over4, prec);
   mpfr_div_2ui (mc->sqrt2_over4, mc->sqrt2_over2, 1ul, GMP_RNDN);

   q24 = (mpc_t *) malloc (mc->cl.h12 * sizeof (mpc_t));
   mc->eta_value = (mpc_t *) malloc (mc->cl.h12 * sizeof (mpc_t));

   for (i = 0; i < mc->cl.h12; i++) {
      mpc_init2 (q24 [i], prec);
      mpc_init2 (mc->eta_value [i], prec);
   }

   if (!checkpoints || !read_q24eta (*mc, q24, "q24")) {
      compute_q24 (*mc, q24, verbose);
      if (checkpoints)
         write_q24eta (*mc, q24, "q24");
   }

   if (!checkpoints || !read_q24eta (*mc, mc->eta_value, "eta")) {
      cm_timer_start (clock2);
      for (i = 0; i < mc->cl.h12; i++) {
         cm_modular_eta_series (mc->m, mc->eta_value [i], q24 [i]);
         if (verbose && i % 200 == 0) {
            printf (".");
            fflush (stdout);
         }
      }
      cm_timer_stop (clock2);
      if (verbose)
         printf ("\n- Time for series:                %.1f\n",
                 cm_timer_get (clock2));
      cm_timer_stop (clock2);
      if (checkpoints)
         write_q24eta (*mc, mc->eta_value, "eta");
   }

   cm_timer_stop (clock1);
   if (verbose)
      printf ("- Time for eta:                   %.1f\n",
              cm_timer_get (clock1));

   for (i = 0; i < mc->cl.h12; i++)
      mpc_clear (q24 [i]);
   free (q24);
}

/*****************************************************************************/

void cm_modclass_clear (cm_modclass_t *mc)

{
   int i;

   mpfr_clear (mc->root);
   for (i = 0; i < mc->cl.h12; i++)
      mpc_clear (mc->eta_value [i]);
   free (mc->eta_value);

   cm_classgroup_clear (&(mc->cl));
   cm_modular_clear (&(mc->m));
}

/*****************************************************************************/
/*                                                                           */
/* file handling functions                                                   */
/*                                                                           */
/*****************************************************************************/

static void write_q24eta (cm_modclass_t mc, mpc_t *q24eta, char *type)
   /* writes the values of q24eta to the file                                */
   /* CLASS_TMPDIR + "/tmp_" + d + "_" + prec + "_" + type + ".dat"          */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%d_%s.dat", CM_CLASS_TMPDIR, -mc.cl.d,
      (int) mc.m.prec, type);

   if (!cm_file_open_write (&f, filename))
      exit (1);

   for (i = 0; i < mc.cl.h12; i++) {
      mpc_out_str (f, 16, 0, q24eta [i], MPC_RNDNN);
      fprintf (f, "\n");
   }

   cm_file_close (f);
}

/*****************************************************************************/

static bool read_q24eta (cm_modclass_t mc, mpc_t *q24eta, char *type)
   /* reads the values of q24eta from the file                               */
   /* CLASS_TMPDIR + "/tmp_" + d + "_" + prec + "_" + type + ".dat"          */

{
   char filename [255];
   FILE *f;
   int i;

   sprintf (filename, "%s/tmp_%"PRIicl"_%d_%s.dat", CM_CLASS_TMPDIR, -mc.cl.d,
            (int) mc.m.prec, type);

   if (!cm_file_open_read (&f, filename))
      return false;

   for (i = 0; i < mc.cl.h12; i++) {
      mpc_inp_str (q24eta [i], f, NULL, 16, MPC_RNDNN);
   }

   cm_file_close (f);

   return true;
}

/*****************************************************************************/
/*                                                                           */
/* functions for precomputations                                             */
/*                                                                           */
/*****************************************************************************/

static void compute_q24 (cm_modclass_t mc, mpc_t *q24, bool verbose)
   /* computes the q^(1/24) for all forms of the class group in mc           */
{
   int i, j, tmp_int;
   mpfr_t Pi24, Pi24_root, tmp;
   mpfr_t *q_real;
   unsigned long int *A_red, *B_red, *order;
   cm_timer  clock2, clock3;
   int    counter1, counter2;

   mpfr_init2 (Pi24, mc.m.prec);
   mpfr_init2 (Pi24_root, mc.m.prec);
   mpfr_init2 (tmp, mc.m.prec);

   mpfr_div_ui (Pi24, mc.m.pi, 24ul, GMP_RNDN);
   mpfr_mul (Pi24_root, Pi24, mc.root, GMP_RNDN);

   A_red = (unsigned long int *) malloc (mc.cl.h12 * sizeof (unsigned long int));
   B_red = (unsigned long int *) malloc (mc.cl.h12 * sizeof (unsigned long int));
   order = (unsigned long int *) malloc (mc.cl.h12 * sizeof (unsigned long int));
   q_real = (mpfr_t *) malloc (mc.cl.h12 * sizeof (mpfr_t));

   for (i = 0; i < mc.cl.h12; i++)
      mpfr_init2 (q_real [i], mc.m.prec);

   cm_timer_start (clock2);
   /* precompute the absolute values of q^{1/24}, i.e. the */
   /* exp (- pi/24 * \sqrt (|d|) / A) for the occurring A  */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   for (i = mc.cl.h12 - 1; i >= 0; i--) {
      /* Check whether the current A is a divisor of a previous one; */
      /* if so, take the previous one which is closest.              */
      /* Notice that the A's are in increasing order.                */
      for (j = i+1;
           j < mc.cl.h12 && mc.cl.form [j].a % mc.cl.form [i].a != 0; j++);
      if (j < mc.cl.h12) {
         if (mc.cl.form [i].a == mc.cl.form [j].a)
            mpfr_set (q_real [i], q_real [j], GMP_RNDN);
         else {
            counter1++;
            mpfr_pow_ui (q_real [i], q_real [j],
                         mc.cl.form [j].a / mc.cl.form [i].a, GMP_RNDN);
         }
      }
      else {
         counter1++;
         counter2++;
         mpfr_div_ui (q_real [i], Pi24_root, mc.cl.form [i].a, GMP_RNDN);
         mpfr_neg (q_real [i], q_real [i], GMP_RNDN);
         mpfr_exp (q_real [i], q_real [i], GMP_RNDN);
      }
      if (verbose &&i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("\n- Number of distinct A:           %d\n", counter1);
      printf ("- Number of exp:                  %d\n", counter2);
      printf ("- Time for exp:                   %.1f\n", cm_timer_get (clock3));
   }

   /* precompute the arguments of q^{1/24}, i.e. the                     */
   /* exp (- i pi/24 * B / A) for the occurring A, B                     */
   /* Write B / A = B_red / A_red with coprime (B_red, A_red).           */
   /* Then order the forms with respect to increasing A_red, and for the */
   /* same A_red with respect to decreasing B_red                        */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   for (i = 0; i < mc.cl.h12; i++) {
      tmp_int = cm_classgroup_gcd (mc.cl.form [i].a, mc.cl.form [i].b);
      A_red [i] = mc.cl.form [i].a / tmp_int;
      B_red [i] = mc.cl.form [i].b / tmp_int;
      order [i] = i;
   }

   /* sort by insertion */
   for (i = 0; i < mc.cl.h12 - 1; i++) {
      tmp_int = i;
      for (j = tmp_int + 1; j < mc.cl.h12; j++)
         if (A_red [order [j]] < A_red [order [tmp_int]]
             || (   A_red [order [j]] == A_red [order [tmp_int]]
             && B_red [order [j]] > B_red [order [tmp_int]]))
            tmp_int = j;
      j = order [i];
      order [i] = order [tmp_int];
      order [tmp_int] = j;
   }

   /* put q24 [] = exp (- pi/24 i / A_red), i.e. the primitive root  */
   /* of unity                                                       */
   for (i = mc.cl.h12 - 1; i >= 0; i--) {
      /* Check whether the current A is a divisor of a previous one; */
      /* if so, take the previous one which is closest.              */
      /* Notice that the A's are in increasing order.                */
      for (j = i+1;
           j < mc.cl.h12 && A_red [order [j]] % A_red [order [i]] != 0;
           j++);
      if (j < mc.cl.h12) {
         if (A_red [order [i]] == A_red [order [j]])
            mpc_set (q24 [order [i]], q24 [order [j]], MPC_RNDNN);
         else {
            counter1++;
            mpc_pow_ui (q24 [order [i]], q24 [order [j]],
                        A_red [order [j]] / A_red [order [i]]);
         }
      }
      else {
         counter1++;
         counter2++;
         mpfr_div_ui (tmp, Pi24, A_red [order [i]], GMP_RNDN);
         mpfr_sin_cos (q24 [order [i]]->im, q24 [order [i]]->re,
                       tmp, GMP_RNDN);
      }
      if (verbose && i % 200 == 0) {
         printf (".");
         fflush (stdout);
      }
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("\n- Number of distinct A_red:       %d\n", counter1);
      printf ("- Number of sin/cos:              %d\n", counter2);
      printf ("- Time for sin/cos:               %.1f\n", cm_timer_get (clock3));
   }

   /* now raise to the power B_red */
   cm_timer_start (clock3);
   counter1 = 0;
   counter2 = 0;
   for (i = mc.cl.h12 - 1; i >= 0; i--)
   {
      if (B_red [order [i]] == 0)
         mpc_set_ui_ui (q24 [order [i]], 1ul, 0ul, MPC_RNDNN);
      else if (B_red [order [i]] > 1) {
         /* Check whether the current B is a multiple of a previous one with */
         /* the same A; if so, take the previous one which is closest.       */
         for (j = i+1;
              j < mc.cl.h12 && A_red [order [j]] == A_red [order [i]]
                    && B_red [order [j]] > 0
                    && B_red [order [i]] % B_red [order [j]] != 0;
              j++);
         if (j < mc.cl.h12 && A_red [order [j]] == A_red [order [i]]
             && B_red [order [j]] > 0) {
            if (B_red [order [i]] == B_red [order [j]])
               mpc_set (q24 [order [i]], q24 [order [j]],
                        MPC_RNDNN);
            else {
               counter1++;
               mpc_pow_ui (q24 [order [i]], q24 [order [j]],
                           B_red [order [i]] / B_red [order [j]]);
            }
         }
         else {
            counter2++;
            mpc_pow_ui (q24 [order[i]], q24 [order [i]],
                        B_red [order [i]]);
         }
      }
   }
   cm_timer_stop (clock3);
   if (verbose) {
      printf ("- Number of partial powers:       %d\n", counter1);
      printf ("- Number of total powers:         %d\n", counter2);
      printf ("- Time for powers:                %.1f\n", cm_timer_get (clock3));
   }

   /* Compute the q^(1/24) in q24 */
   for (i = 0; i < mc.cl.h12; i++)
      mpc_mul_fr (q24 [i], q24 [i], q_real [i], MPC_RNDNN);
   cm_timer_stop (clock2);
   if (verbose)
      printf ("- Time for q^(1/24):              %.1f\n", cm_timer_get (clock2));

   for (i = 0; i < mc.cl.h12; i++)
      mpfr_clear (q_real [i]);

   free (A_red);
   free (B_red);
   free (order);
   free (q_real);

   mpfr_clear (Pi24);
   mpfr_clear (Pi24_root);
   mpfr_clear (tmp);
}

/*****************************************************************************/
/*                                                                           */
/* other internal functions                                                  */
/*                                                                           */
/*****************************************************************************/

static void cm_modclass_mpc_set_quadratic (cm_modclass_t mc, mpc_t rop,
                                        int_cl_t a, int_cl_t b)
      /* sets rop to (b + sqrt (d)) / (2a)                                      */

{
   mpfr_set_si (rop->re, b, GMP_RNDN);
   mpfr_set (rop->im, mc.root, GMP_RNDN);
   mpc_div_ui (rop, rop, 2*a, MPC_RNDNN);
}

/*****************************************************************************/

static void cm_modclass_fundamental_domain_quad (int_cl_t d, int_cl_t *a,
   int_cl_t *b, cm_matrix_t *M)
      /* transforms (b + sqrt (d)) / (2a) into the fundamental domain and   */
      /* returns the inverse transformation matrix M                        */
      /* Unfortunately we have to compute internally with multiprecision;   */
      /* although the final result fits into a long, intermediate results   */
      /* can be quite large.                                                */
{
   bool reduced = false;
   static mpz_t a_local, b_local, b_minus_a, two_a, offset, c;
   static bool first_time = true;
   int_cl_t     tmp_int;

   if (first_time)
   {
      mpz_init (a_local);
      mpz_init (b_local);
      mpz_init (b_minus_a);
      mpz_init (two_a);
      mpz_init (offset);
      mpz_init (c);
      first_time = false;
   }

   M->a = 1;
   M->b = 0;
   M->c = 0;
   M->d = 1;
   mpz_set_si (a_local, *a);
   mpz_set_si (b_local, *b);

   while (!reduced) {
      /* obtain -a < b <= a                                                   */
      /* basically, compute offset as b / (2a) rounded towards the closest    */
      /* integer. If the quotient is exactly between two integers, round down */
      mpz_sub (b_minus_a, b_local, a_local);
      mpz_mul_2exp (two_a, a_local, 1);
      mpz_cdiv_q (offset, b_minus_a, two_a);

      tmp_int = mpz_get_si (offset);
      /* multiply M from the right by T^{tmp_int} */
      M->b += M->a * tmp_int;
      M->d += M->c * tmp_int;

      mpz_mul (offset, offset, two_a);
      mpz_sub (b_local, b_local, offset);

      /* compute c */
      mpz_mul (c, b_local, b_local);
      mpz_sub_si (c, c, d);
      mpz_fdiv_q_2exp (c, c, 2);
      mpz_divexact (c, c, a_local);
      /* if not reduced, invert */
      if (mpz_cmp (a_local, c) < 0 ||
          (mpz_cmp (a_local, c) == 0 && mpz_sgn (b_local) > 0))
         reduced = true;
      else {
         mpz_set (a_local, c);
         mpz_neg (b_local, b_local);
         tmp_int = M->a;
         M->a = M->b;
         M->b = -tmp_int;
         tmp_int = M->c;
         M->c = M->d;
         M->d = -tmp_int;
         reduced = false;
      }
   }

   *a = mpz_get_si (a_local);
   *b = mpz_get_si (b_local);
   /* normalise the matrix */
   if (M->c < 0 || (M->c == 0 && M->d < 0)) {
      M->a = -M->a;
      M->b = -M->b;
      M->c = -M->c;
      M->d = -M->d;
   }
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_modclass_eta_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b)
   /* evaluates the eta function at the quadratic integer                    */
   /* (b + sqrt (d)) / (2a)                                                  */

{
   int_cl_t a_local, b_local;
   int      i, sign;
   cm_matrix_t M;
   mpc_t    tmp;

   a_local = a;
   b_local = b;
   cm_modclass_fundamental_domain_quad (mc.cl.d, &a_local, &b_local, &M);
   cm_modclass_mpc_set_quadratic (mc, rop, a_local, b_local);
   cm_modular_eta_transform (mc.m, rop, rop, M);

   /* look up the eta value */
   i = 0;
   if (b_local < 0)
   {
      b_local = -b_local;
      sign = -1;
   }
   else
      sign = 1;
   while ((i < mc.cl.h12)
     && (mc.cl.form [i].a != a_local || mc.cl.form [i].b != b_local))
      i++;
   if (i == mc.cl.h12)
   {
      /* eta value not found, compute it. May happen when the level of the   */
      /* modular function and the conductor have a common factor. This case  */
      /* is rare, and the following computations are not optimised.          */
      printf ("Q");
      mpc_init2 (tmp, mpc_get_prec (rop));
      if (sign == 1)
         cm_modclass_mpc_set_quadratic (mc, tmp, a_local, b_local);
      else
         cm_modclass_mpc_set_quadratic (mc, tmp, a_local, -b_local);
      cm_modular_eta_eval (mc.m, tmp, tmp);
      mpc_mul (rop, rop, tmp, MPC_RNDNN);
      mpc_clear (tmp);
   }
   else
   {
      if (sign == 1)
         mpc_mul (rop, rop, mc.eta_value [i], MPC_RNDNN);
      else
      {
         mpc_conj (mc.eta_value [i], mc.eta_value [i], MPC_RNDNN);
         mpc_mul (rop, rop, mc.eta_value [i], MPC_RNDNN);
         mpc_conj (mc.eta_value [i], mc.eta_value [i], MPC_RNDNN);
      }
   }
}

/*****************************************************************************/

void cm_modclass_f_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b)
   /* evaluates the Weber f-function at the quadratic integer                */
   /* (b + sqrt (d)) / (2*a)                                                 */
   /* We get into trouble since the numerator sometimes belongs to a         */
   /* non-primitive form, i.e. a form corresponding to d/4.                  */
   /* Even worse, things become completely mixed up when d/4 is not a        */
   /* discriminant; there, one might have to use forms with half-integral    */
   /* a and b.                                                               */
   /* We obtain the denominator from the precomputed eta values and          */
   /* recompute the numerator.                                               */

{
   mpc_t z, tmp;

   mpc_init2 (z, mpc_get_prec (rop));
   mpc_init2 (tmp, mpc_get_prec (rop));

   cm_modclass_eta_eval_quad (mc, tmp, a, b);
   cm_modclass_mpc_set_quadratic (mc, z, a, b);
   mpc_add_ui (z, z, 1ul, MPC_RNDNN);
   mpc_div_ui (z, z, 2ul, MPC_RNDNN);
   cm_modular_eta_eval (mc.m, rop, z);
   mpc_div (rop, rop, tmp, MPC_RNDNN);
   mpc_mul (rop, rop, mc.m.zeta48inv, MPC_RNDNN);

   mpc_clear (z);
   mpc_clear (tmp);
}

/*****************************************************************************/

void cm_modclass_f1_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b)
   /* evaluates the Weber f1-function at the quadratic integer               */
   /* (b + sqrt (d)) / (2*a); the same comments as for f apply.              */

{
   mpc_t z, tmp;

   mpc_init2 (z, mpc_get_prec (rop));
   mpc_init2 (tmp, mpc_get_prec (rop));

   cm_modclass_eta_eval_quad (mc, rop, a, b);
   cm_modclass_mpc_set_quadratic (mc, z, a, b);
   mpc_div_ui (tmp, z, 2ul, MPC_RNDNN);
   cm_modular_eta_eval (mc.m, tmp, tmp);
   mpc_div (rop, tmp, rop, MPC_RNDNN);

   mpc_clear (z);
   mpc_clear (tmp);
}

/*****************************************************************************/

void cm_modclass_j_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b)

{
   mpc_t f;

   mpc_init2 (f, mpc_get_prec (rop));

   cm_modclass_f1_eval_quad (mc, f, a, b);
   mpc_pow_ui (f, f, 8ul);
   mpc_sqr (rop, f, MPC_RNDNN);
   mpc_ui_div (f, 16ul, f, MPC_RNDNN);
   mpc_add (rop, rop, f, MPC_RNDNN);
   mpc_pow_ui (rop, rop, 3ul);

   mpc_clear (f);
}

/*****************************************************************************/

void cm_modclass_gamma2_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b)

{
   mpc_t f;

   mpc_init2 (f, mpc_get_prec (rop));

   cm_modclass_f1_eval_quad (mc, f, a, b);
   mpc_pow_ui (f, f, 8ul);
   mpc_sqr (rop, f, MPC_RNDNN);
   mpc_ui_div (f, 16ul, f, MPC_RNDNN);
   mpc_add (rop, rop, f, MPC_RNDNN);

   mpc_clear (f);
}

/*****************************************************************************/

void cm_modclass_gamma3_eval_quad (cm_modclass_t mc, mpc_t rop,
   int_cl_t a, int_cl_t b)
   /* evaluates sqrt (d) * gamma3 */

{
   mpc_t f, tmp;
   mpfr_t tmp_fr;

   mpc_init2 (f, mpc_get_prec (rop));
   mpc_init2 (tmp, mpc_get_prec (rop));

   cm_modclass_f_eval_quad (mc, f, a, b);
   mpc_pow_ui (f, f, 8ul);
   cm_modclass_f1_eval_quad (mc, rop, a, b);
   mpc_pow_ui (rop, rop, 8ul);

   mpc_mul_ui (rop, rop, 2ul, MPC_RNDNN);
   mpc_sub (rop, rop, f, MPC_RNDNN);
   mpc_pow_ui (tmp, f, 3ul);
   mpc_add_ui (tmp, tmp, 8ul, MPC_RNDNN);
   mpc_mul (rop, rop, tmp, MPC_RNDNN);
   mpc_div (rop, rop, f, MPC_RNDNN);
   mpc_mul_fr (rop, rop, mc.root, MPC_RNDNN);
   /* multiply by i */
   tmp_fr [0] = rop->im [0];
   rop->im [0] = rop->re [0];
   rop->re [0] = tmp_fr [0];
   mpfr_neg (rop->re, rop->re, GMP_RNDN);

   mpc_clear (f);
   mpc_clear (tmp);
}

/*****************************************************************************/
