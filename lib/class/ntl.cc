#include "cm_class-impl.h"

#include <sstream>
#include <NTL/ZZ_pXFactoring.h>

using namespace NTL;
using namespace std;

/*****************************************************************************/

static void mpz_get_ZZ (ZZ &result, mpz_t arg)

{
   stringstream tmp;
   char *       tmp2;

   tmp2 = mpz_get_str (NULL, 10, arg);

   tmp << tmp2;
   tmp >> result;

   free (tmp2);
}

/*****************************************************************************/

static void ZZp_get_mpz (mpz_t result, ZZ_p arg)

{
   ostringstream tmp;

   tmp << arg << ends;
   mpz_set_str (result, tmp.str ().c_str (), 10);
}

/*****************************************************************************/

static void mpz_get_ZZp (ZZ_p &result, mpz_t arg)

{
   stringstream tmp;
   char *       tmp2;

   tmp2 = mpz_get_str (NULL, 10, arg);

   tmp << tmp2;
   tmp >> result;

   free (tmp2);
}

/*****************************************************************************/

void cm_ntl_find_factor (mpz_t *res, mpz_t *f_z, int f_deg, int factor_deg,
   mpz_t p, bool verbose)
   // assumes that the polynomial f is monic of degree f_deg and splits into
   // distinct factors of degree factor_deg over the prime field
   // One of the factors is returned as a list of coefficients via res, which
   // must contain sufficiently many initialised fields.
   // We use a trick by Atkin to speed up the splitting process. If n | p-1,
   // then there is a primitive n-th root of unity zeta in the prime field.
   // We may compute the gcd of all the X^((p^degree-1) / n) - zeta^i with f
   // to break f into smaller pieces. To do so, only one exponentiation modulo
   // f is needed.

{
   mpz_t             exponent, exponent_prime, primitive, zeta;
   mpz_t             tmp_z;
   ZZ                p_ZZ, tmp_ZZ, exponent_ZZ;
   mpz_get_ZZ (p_ZZ, p);
   ZZ_p::init (p_ZZ);
   ZZ_p              zeta_ZZ;
   unsigned long int n, n_factors [17];
   uint_cl_t         n_factors_cl [17];
   unsigned int      n_exp [17];
      // n and its prime factorisation
   bool              ok;
   int               i, start_degree, min_degree, min_degree_arg;
   ZZ_p              a, zeta_power;
   ZZ_pX             *factors;
   ZZ_pX             powerx, powerx_zeta;
   ZZ_pXModulus      F;
   cm_timer             clock, clock2;

   cm_timer_start (clock2);
   if (verbose)
      cout << "--- Root finding, degree " << factor_deg << " of " << f_deg << endl;
   mpz_init (exponent);
   mpz_init (exponent_prime);
   mpz_init (primitive);
   mpz_init (zeta);
   mpz_init (tmp_z);

   mpz_pow_ui (exponent, p, factor_deg);
   mpz_sub_ui (exponent, exponent, 1);
   mpz_sub_ui (exponent_prime, p, 1);

   // We are looking for an n-th root of unity in the prime field, even when
   // factors are sought over an extension field.
   // look for an n slightly smaller than twice the degree
   for (n = 2 * f_deg; mpz_divisible_ui_p (exponent_prime, n) == 0; n--) ;
   mpz_divexact_ui (exponent_prime, exponent_prime, n);
   cm_classgroup_factor ((int_cl_t) n, n_factors_cl, n_exp);
   for (i = 0; i < 17; i++)
      n_factors [i] = (unsigned long int) n_factors_cl [i];

   mpz_divexact_ui (exponent, exponent, n);
   if (verbose)
      cout << "n " << n << endl;

   // compute a primitive n-th root of unity in the prime field
   if (n == 1)
      mpz_set_ui (zeta, 1);
   else {
      mpz_set_ui (primitive, 1);
      ok = false;
      while (!ok) {
         mpz_add_ui (primitive, primitive, 1);
         mpz_powm (zeta, primitive, exponent_prime, p);
         if (mpz_cmp_ui (zeta, 1) != 0) {
            ok = true;
            i = 0;
            while (ok && n_factors [i] != 0) {
               if (n_exp [i] != 0) {
                  mpz_powm_ui (tmp_z, zeta, n / n_factors [i], p);
                  if (mpz_cmp_ui (tmp_z, 1) == 0)
                     ok = false;
               }
               i++;
            }
         }
      }
   }

   // copy the polynomial into a variable of NTL type
   ZZ_pX f;
   f.rep.SetLength (f_deg + 1);
   for (i = 0; i <= f_deg; i++)
      mpz_get_ZZp (f.rep [i], f_z [i]);

   mpz_get_ZZ (exponent_ZZ, exponent);
   mpz_get_ZZp (zeta_ZZ, zeta);

   factors = new ZZ_pX [n];
   factors [0] = f;
   a = 0;
   while (deg (factors [0]) > factor_deg) {
      // Split the polynomial in factors [0]. Notice that after each execution
      // of the code in the loop, the product of factors [0], ..., factors [i]
      // remains the same. At the same time, we keep track of the factor of
      // minimal degree (not including factor [0]).
      start_degree = deg (factors [0]);
      if (verbose)
         cout << "--- degree " << start_degree << endl;
      cm_timer_start (clock);
      build (F, factors [0]);
      min_degree = start_degree + 1;
      min_degree_arg = -1;
      PowerXPlusAMod (powerx, a, exponent_ZZ, F);
      cm_timer_stop (clock);
      if (verbose)
         cout << "-- Time for power " << cm_timer_get (clock) << endl;
      zeta_power = 1;
      cm_timer_start (clock);
      for (i = 1; (unsigned int) i < n && deg (factors [0]) > factor_deg
                  && min_degree > (1.5 * start_degree) / n
                  && min_degree > factor_deg; i++)
      // Stop as soon as one of the factors has the desired degree or
      // factors [0] = 1.
      // The expected outcome of one gcd step is a factor of degree (original
      // degree divided by n). To avoid computing too many gcds, we stop as
      // soon as we are close to this value
      {
         mul (zeta_power, zeta_power, zeta_ZZ);
         sub (powerx_zeta, powerx, zeta_power);
         GCD (factors [i], powerx_zeta, factors [0]);
         div (factors [0], factors [0], factors [i]);
         if (deg (factors [i]) > 0 && deg (factors [i]) < min_degree)
         {
            min_degree = deg (factors [i]);
            min_degree_arg = i;
         }
         if (verbose) {
            if (deg (factors [i]) > 0)
               cout << deg (factors [i]) << " " << flush;
            else
               cout << ". " << flush;
         }
      }
      // In two cases, the splitting was unsuccesful: If all the factors [i],
      // i >= 1, are 1 (which means min_degree_arg equals -1) or if one of
      // them contains the complete original factors [0], which means that
      // min_degree equals start_degree.
      cm_timer_stop (clock);
      if (verbose)
         cout << endl << "-- Time for gcd " << cm_timer_get (clock) << endl;
      if (min_degree_arg == -1) {
         if (verbose)
            cout << a << flush;
      }
      else if (min_degree == start_degree) {
         factors [0] = factors [min_degree_arg];
         if (verbose)
            cout << a << flush;
      }
      else {
         // The real factor of minimal degree might be in factors [0],
         // otherwise swap.
         if (deg (factors [0]) == 0 || deg (factors [0]) > min_degree)
            factors [0] = factors [min_degree_arg];
         if (verbose)
            cout << "." << flush;
      }
      a++;
   }

   if (verbose)
      cout << endl;
   // copy the result from factors [0]
   for (i = 0; i <= factor_deg; i++)
      ZZp_get_mpz (res [i], factors [0].rep [i]);
   delete [] factors;
   mpz_clear (exponent);
   mpz_clear (exponent_prime);
   mpz_clear (primitive);
   mpz_clear (zeta);
   mpz_clear (tmp_z);
   cm_timer_stop (clock2);
   if (verbose)
      cout << "-- Time for root " << cm_timer_get (clock2) << endl;
}

/*****************************************************************************/

void cm_ntl_find_root_split (mpz_t root, mpz_t *f, int deg, mpz_t p,
   bool verbose)
   // finds a root of the monic polynomial f of degree deg over the prime
   // field of characteristic p
   // assumes that the polynomial splits into distinct linear factors
   // root is changed

{
   mpz_t factor [2];

   mpz_init (factor [0]);
   mpz_init (factor [1]);
   cm_ntl_find_factor (factor, f, deg, 1, p, verbose);
   if (mpz_sgn (factor [0]) != 0)
      mpz_sub (root, p, factor [0]);
   else
      mpz_set_ui (root, 0);
   mpz_clear (factor [0]);
   mpz_clear (factor [1]);
}

/*****************************************************************************/

mpz_t* cm_ntl_find_roots (mpz_t *f_z, int f_deg, mpz_t p, int *no)
   // computes all the roots of the polynomial f of degree deg in the prime
   // field of characteristic p without multiplicities.
   // no is the number of found roots
   // Since this will usually be called with a polynomial of small degree,
   // the implementation does not use any tricks.

{
   ZZ           exponent;
   int          i, no_old, no_new;
   ZZ           p_ZZ;
   mpz_get_ZZ (p_ZZ, p);
   ZZ_p::init (p_ZZ);
   ZZ_p         a;
   mpz_t        *result;
   ZZ_pX        f, x (1, 1), powerx, *factors;
   ZZ_pXModulus F;
   cm_timer           clock, clock2;

   cm_timer_start (clock);
//   cout << "--- Root finding" << endl;
   // copy the polynomial into a variable of NTL type
   f.rep.SetLength (f_deg + 1);
   for (i = 0; i <= f_deg; i++)
      mpz_get_ZZp (f.rep [i], f_z [i]);

   factors = new ZZ_pX [f_deg];

   // isolate the part with roots in the prime field
   MakeMonic (f);
   build (F, f);
   cm_timer_start (clock2);
   PowerXMod (powerx, p_ZZ, F);
   cm_timer_stop (clock2);
//   cout << "--- Time for power " << cm_timer_get (clock2) << endl;
   sub (powerx, powerx, x);
   GCD (factors [0], powerx, f);
   *no = deg (factors [0]);
   a = 0;
   if (*no > 0)
   {
      build (F, factors [0]);
      sub (exponent, p_ZZ, 1);
      div (exponent, exponent, 2);
      no_old = 1;
      no_new = 1;
         // the number of found factors at the beginning and at the end
         // of the loop
      while (no_old < *no)
      {
         for (i = 0; i < no_old; i++)
            if (deg (factors [i]) > 1)
            // split the polynomial in factors [i]
            {
               PowerXPlusAMod (powerx, a, exponent, factors [i]);
               sub (powerx, powerx, 1);
               GCD (factors [no_new], powerx, factors [i]);
               if (   deg (factors [no_new]) > 0
                   && deg (factors [no_new]) < deg (factors [i]))
               {
                  div (factors [i], factors [i], factors [no_new]);
                  no_new++;
               }
            }
         no_old = no_new;
         a++;
      }
   }

   result =(mpz_t*) malloc ((*no) * sizeof (mpz_t));
   for (i = 0; i < *no; i++)
   {
      mpz_init (result [i]);
      ZZp_get_mpz (result [i], factors [i].rep [0]);
      if (mpz_sgn (result [i]) != 0)
         mpz_sub (result [i], p, result [i]);
   }

   cm_timer_stop (clock);
   delete [] factors;

   return result;
}

/*****************************************************************************/
