#include "cm_common-impl.h"

#define MAX(a,b) ((a > b ? a : b))

static bool find_in_chain (int* index, cm_qdev_t f, int length, long int no);
static long int lognormmax (mpc_t op);
static int qdev_atkin_coeff (int i, int l);
static int qdev_atkinhecke_coeff (int i, int l, int r);

/*****************************************************************************/
/*                                                                           */
/* internal functions                                                        */
/*                                                                           */
/*****************************************************************************/

static bool find_in_chain (int* index, cm_qdev_t f, int length, long int no)
   /* looks for no in the first length elements of f                         */
   /* The return value indicates the success; if the operation succeeds,     */
   /* the index contains i such that f.chain [i][0] = no.                    */

{
   int left = 0, right = length - 1, middle;

   if (no < f.chain [0][0] || no > f.chain [length-1][0])
      return false;

   while (left < right - 1)
   {
      middle = (left + right) / 2;
      if (f.chain [middle][0] < no)
         left = middle;
      else
         right = middle;
   }
   if (f.chain [left][0] == no)
   {
      *index = left;
      return true;
   }
   else if (f.chain [right][0] == no)
   {
      *index = right;
      return true;
   }
   else
      return false;
}

/*****************************************************************************/

static long int lognormmax (mpc_t op)
   /* computes the logarithm in base 2 (as the exponent of the value) of the */
   /* max norm of op                                                         */

{
   if (mpfr_sgn (op->re) == 0)
      if (mpfr_sgn (op->im) == 0)
         return mpfr_get_emin ();
      else
         return (op->im->_mpfr_exp);
   else
      if (mpfr_sgn (op->im) == 0)
         return (op->re->_mpfr_exp);
      else
         return MAX(op->re->_mpfr_exp, op->im->_mpfr_exp);
}

/*****************************************************************************/

static int qdev_atkin_coeff (int i, int l)
   /* returns the coefficient in front of (q^(1/24))^i of the Atkin function */
   /* eta (z) * (eta (l z)                                                   */

{
   long int n;
   long int m;
   int      k, c;

   if (i < l+1 || (i - (l+1)) % 24 != 0)
      return 0;

   c = 0;
   k = (i - (l+1)) / 24;

   for (n = 0; l*n*(3*n-1) <= 2*k; n++)
   {
      /* Let m be the integral solution of m (3m-1)/2 = k - l n (3n-1) / 2,  */
      /* if it exists.                                                       */
      m = cm_nt_sqrt (1 + 24*k - 12*l*n*(3*n-1));
      if (m != -1)
      {
         if ((m - 1) % 6 == 0)
            m = (1 - m) / 6;
         else
            m = (1 + m) / 6;
         if ((m + n) % 2 == 0)
            c += 1;
         else
            c -= 1;
      }
   }
   for (n = 1; l*n*(3*n+1) <= 2*k; n++)
   {
      m = cm_nt_sqrt (1 + 24*k - 12*l*n*(3*n+1));
      if (m != -1)
      {
         if ((m - 1) % 6 == 0)
            m = (1 - m) / 6;
         else
            m = (1 + m) / 6;
         if ((m + n) % 2 == 0)
            c += 1;
         else
            c -= 1;
      }
   }

   return c;
}

/*****************************************************************************/

static int qdev_atkinhecke_coeff (int i, int l, int r)
   /* Let the r-th Hecke operator applied to the Atkin function evaluated in */
   /* q^24 (that is, eta (24 z) * eta (24 l z)) result in a q-series         */
   /* q^(l+1) \sum c_i q^(24 i). Then T_r applied to eta (z) * eta (lz)      */
   /* is q^((l+1)/24) \sum c_i q^i.                                          */
   /* The function returns c_i.                                              */


{
   int I = 24*i + l+1;
   if (I % r != 0)
      return qdev_atkin_coeff (I*r, l);
   else
      return qdev_atkin_coeff (I*r, l) + qdev_atkin_coeff (I/r, l);
}

/*****************************************************************************/
/*                                                                           */
/* external functions                                                        */
/*                                                                           */
/*****************************************************************************/

void cm_qdev_init (cm_qdev_t *f)
   /* initialises the addition chain for eta */

{
   int n, i, j, k;

   f->length = 321;
   /* must be odd                                                        */
   /* 321 corresponds to about 300000 bits in the worst case, each power */
   /* of q yielding at least log_2 (exp (sqrt (3) * pi)) = 7.85 bits.    */

   f->chain = (long int **) malloc (f->length * sizeof (long int *));
   for (n = 0; n < f->length; n++)
      f->chain [n] = (long int *) malloc (6 * sizeof (long int));

   (*f).chain [0][0] = 0;
   (*f).chain [0][5] = 1;
   for (n = 1; n <= f->length / 2; n++)
   {
      (*f).chain [2*n-1][0] = n*(3*n-1) / 2;
      (*f).chain [2*n][0] = n*(3*n+1) / 2;
      if (n % 2 == 0)
      {
         (*f).chain [2*n-1][5] = 1;
         (*f).chain [2*n][5] = 1;
      }
      else
      {
         (*f).chain [2*n-1][5] = -1;
         (*f).chain [2*n][5] = -1;
      }
   }

   (*f).chain [0][1] = 0;
   (*f).chain [1][1] = 0;
   for (n = 2; n < f->length; n++)
   {
      (*f).chain [n][1] = 0;
      /* try to express an even exponent as twice a previous one      */
      if ((*f).chain [n][0] % 2 == 0)
         if (find_in_chain (&i, *f, n, (*f).chain [n][0] / 2))
            {
               (*f).chain [n][1] = 1;
               (*f).chain [n][2] = i;
            }
      /* try to express the exponent as the sum of two previous ones */
      for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
         if (find_in_chain (&j, *f, n, (*f).chain [n][0] - (*f).chain [i][0]))
            {
               (*f).chain [n][1] = 2;
               (*f).chain [n][2] = i;
               (*f).chain [n][3] = j;
            }
      /* try to express an exponent as four times a previous one      */
      if ((*f).chain [n][0] % 4 == 0)
         if (find_in_chain (&i, *f, n, (*f).chain [n][0] / 4))
            {
               (*f).chain [n][1] = 3;
               (*f).chain [n][2] = i;
            }
      /* try to express an even exponent as twice the sum of two previous */
      /* ones */
      /* try to express the exponent as twice a previous plus a third one */
      for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
         if (find_in_chain (&j, *f, n, (*f).chain [n][0] - 2 * (*f).chain [i][0]))
      {
         (*f).chain [n][1] = 4;
         (*f).chain [n][2] = i;
         (*f).chain [n][3] = j;
      }
      if ((*f).chain [n][0] % 2 == 0)
         for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
            if (find_in_chain (&j, *f, n, (*f).chain [n][0]/2 - (*f).chain [i][0]))
               {
                  (*f).chain [n][1] = 5;
                  (*f).chain [n][2] = i;
                  (*f).chain [n][3] = j;
               }
      /* try to express the exponent as the sum of three previous ones */
      for (i = 0; i < n && (*f).chain [n][1] == 0; i++)
         for (j = i; j < n && (*f).chain [n][1] == 0; j++)
            if (find_in_chain (&k, *f, n, (*f).chain [n][0] - (*f).chain [i][0]
                 - (*f).chain [j][0]))
            {
               (*f).chain [n][1] = 6;
               (*f).chain [n][2] = i;
               (*f).chain [n][3] = j;
               (*f).chain [n][4] = k;
            }
      if ((*f).chain [n][1] == 0)
      {
         printf ("*** Houston, we have a problem! No success for element ");
         printf ("%i = %li in the addition chain ", n, (*f).chain [n][0]);
         printf ("computation in qdev_init. ");
         printf ("Go back programming!\n");
         exit (1);
      }
   }
}

/*****************************************************************************/

void cm_qdev_clear (cm_qdev_t *f)

{
   int n;

   for (n = 0; n < f->length; n++)
      free (f->chain [n]);
   free (f->chain);
}

/*****************************************************************************/

void cm_qdev_eval (mpc_t rop, cm_qdev_t f, mpc_t q1)
   /* evaluates f in q */

{
   mp_exp_t prec;
   mpc_t    *q, term;
   int      n, i;

   prec = mpfr_get_prec (rop->re);

   q = (mpc_t *) malloc (f.length * sizeof (mpc_t));
   mpc_init2 (q [1], prec);
   mpc_set (q [1], q1, MPC_RNDNN);
   mpc_init2 (term, prec);

   mpc_set_si (rop, f.chain [0][5], MPC_RNDNN);
   if (f.chain [1][5] != 0)
   {
      mpc_mul_si (term, q [1], f.chain [1][5], MPC_RNDNN);
      mpc_add (rop, rop, term, MPC_RNDNN);
   }
   n = 1;

   /* take next power of q into account if result is not precise enough */
   while (lognormmax (q [n]) > -prec)
   {
      n++;
      if (n >= f.length)
      {
         printf ("*** Houston, we have a problem! Addition chain too short ");
         printf ("in 'qdev_eval'.\n");
         printf ("n=%i, length=%i\n", n, f.length);
         printf ("q "); mpc_out_str (stdout, 10, 10, q [1], MPC_RNDNN);
         printf ("\n");
         printf ("q^i "); mpc_out_str (stdout, 10, 10, q [n-1], MPC_RNDNN);
         printf ("\n");
         exit (1);
      }
      mpc_init2 (q [n], prec);
      switch (f.chain [n][1])
      {
      case 1:
         mpc_sqr (q [n], q [f.chain [n][2]], MPC_RNDNN);
         break;
      case 2:
         mpc_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  MPC_RNDNN);
         break;
      case 3:
         mpc_sqr (q [n], q [f.chain [n][2]], MPC_RNDNN);
         mpc_sqr (q [n], q [n], MPC_RNDNN);
         break;
      case 4:
         mpc_sqr (q [n], q [f.chain [n][2]], MPC_RNDNN);
         mpc_mul (q [n], q[n], q [f.chain [n][3]], MPC_RNDNN);
         break;
      case 5:
         mpc_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  MPC_RNDNN);
         mpc_sqr (q [n], q[n], MPC_RNDNN);
         break;
      case 6:
         mpc_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  MPC_RNDNN);
         mpc_mul (q [n], q [n], q [f.chain [n][4]], MPC_RNDNN);
         break;
      }
      if (f.chain [n][5] != 0)
      {
         mpc_mul_si (term, q [n], f.chain [n][5], MPC_RNDNN);
         mpc_add (rop, rop, term, MPC_RNDNN);
      }
   }

   for (i = 1; i <= n; i++)
      mpc_clear (q [i]);
   free (q);
   mpc_clear (term);
}

/*****************************************************************************/

void cm_qdev_eval_fr (mpfr_t rop, cm_qdev_t f, mpfr_t q1)
   /* evaluates f in q */

{
   mp_exp_t prec;
   mpfr_t    *q, term;
   int      n, i;

   prec = mpfr_get_prec (rop);

   q = (mpfr_t *) malloc (f.length * sizeof (mpfr_t));
   mpfr_init2 (q [1], prec);
   mpfr_set (q [1], q1, GMP_RNDN);
   mpfr_init2 (term, prec);

   mpfr_set_si (rop, f.chain [0][5], GMP_RNDN);
   mpfr_mul_si (term, q [1], f.chain [1][5], GMP_RNDN);
   mpfr_add (rop, rop, term, GMP_RNDN);
   n = 1;

   /* take next power of q into account if result is not precise enough */
   while (mpfr_get_exp (term) > -prec)
   {
      n++;
      if (n >= f.length)
      {
         printf ("*** Houston, we have a problem! Addition chain too short ");
         printf ("in 'qdev_eval_fr'.\n");
         printf ("n=%i, length=%i\n", n, f.length);
         printf ("q "); mpfr_out_str (stdout, 10, 0, q [1], GMP_RNDN);
         printf ("\n");
         printf ("q^i "); mpfr_out_str (stdout, 10, 0, q [n-1], GMP_RNDN);
         printf ("\n");
         exit (1);
      }
      mpfr_init2 (q [n], prec);
      switch (f.chain [n][1])
      {
      case 1:
         mpfr_sqr (q [n], q [f.chain [n][2]], GMP_RNDN);
         break;
      case 2:
         mpfr_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  GMP_RNDN);
         break;
      case 3:
         mpfr_sqr (q [n], q [f.chain [n][2]], GMP_RNDN);
         mpfr_sqr (q [n], q [n], GMP_RNDN);
         break;
      case 4:
         mpfr_sqr (q [n], q [f.chain [n][2]], GMP_RNDN);
         mpfr_mul (q [n], q[n], q [f.chain [n][3]], GMP_RNDN);
         break;
      case 5:
         mpfr_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  GMP_RNDN);
         mpfr_sqr (q [n], q[n], GMP_RNDN);
         break;
      case 6:
         mpfr_mul (q [n], q [f.chain [n][2]], q [f.chain [n][3]],
                  GMP_RNDN);
         mpfr_mul (q [n], q [n], q [f.chain [n][4]], GMP_RNDN);
         break;
      }
      if (f.chain [n][5] != 0)
      {
         mpfr_mul_si (term, q [n], f.chain [n][5], GMP_RNDN);
         mpfr_add (rop, rop, term, GMP_RNDN);
      }
   }

   for (i = 1; i <= n; i++)
      mpfr_clear (q [i]);
   free (q);
   mpfr_clear (term);
}

/*****************************************************************************/
/*****************************************************************************/

int cm_qdev_atkinhecke_order (int l, int r)
   /* In the notation of the previous function, returns the first i such     */
   /* c_i is different from 0. This is also the order at infinity of         */
   /* T_r (f) / f for f (z) = eta (z) * eta (l z).                           */
   /* Returns 0 if T_r (f) / f = 1.                                          */

{
   int i;
   for (i = (-(l+1)) / 24; i < 0 && qdev_atkinhecke_coeff (i, l, r) == 0; i++);

   return i;
}

/*****************************************************************************/

