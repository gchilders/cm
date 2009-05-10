#include "mpfp-impl.h"

/*****************************************************************************/

static mpz_t        __mpfp_p;
static unsigned int __mpfp_mul, __mpfp_sqr, __mpfp_inv;
   /* count the number of multiplications, squarings and inversions in the */
   /* field */

/*****************************************************************************/
/*****************************************************************************/

void mpfp_type_init (mpz_t p)

{
   mpz_init_set (__mpfp_p, p);
   __mpfp_mul = 0;
   __mpfp_sqr = 0;
   __mpfp_inv = 0;
}

/*****************************************************************************/
/*****************************************************************************/

void mpfp_init (mpfp_t rop)

{
   mpz_init (rop);
}

/*****************************************************************************/

void mpfp_clear (mpfp_t rop)

{
   mpz_clear (rop);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfp_out (mpfp_t op)

{
   mpz_out_str (stdout, 10, op);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfp_set (mpfp_t rop, mpfp_t op)

{
   mpz_set (rop, op);
}

/*****************************************************************************/

void mpfp_set_z (mpfp_t rop, mpz_t op)

{
   mpz_mod (rop, op, __mpfp_p);
}

/*****************************************************************************/

void mpfp_set_ui (mpfp_t rop, long unsigned int op)

{
   mpz_set_ui (rop, op);
   mpz_mod (rop, rop, __mpfp_p);
}

/*****************************************************************************/

void mpfp_set_str (mpfp_t rop, char *op, int base)

{
   mpz_set_str (rop, op, base);
   mpz_mod (rop, rop, __mpfp_p);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfp_init_set (mpfp_t rop, mpfp_t op)

{
   mpfp_init (rop);
   mpfp_set (rop, op);
}

/*****************************************************************************/

void mpfp_init_set_z (mpfp_t rop, mpz_t op)

{
   mpfp_init (rop);
   mpfp_set_z (rop, op);
}

/*****************************************************************************/

void mpfp_init_set_ui (mpfp_t rop, unsigned long int op)

{
   mpfp_init (rop);
   mpfp_set_ui (rop, op);
}

/*****************************************************************************/
/*****************************************************************************/

void mpfp_get_z (mpz_t rop, mpfp_t op)

{
   mpz_set (rop, op);
}

/*****************************************************************************/
/*****************************************************************************/

int mpfp_cmp_ui (mpfp_t op1, unsigned long int op2)

{
   return mpz_cmp_ui (op1, op2);
}

/*****************************************************************************/

int mpfp_cmp_si (mpfp_t op1, long int op2)

{
   return mpz_cmp_si (op1, op2);
}

/*****************************************************************************/

void mpfp_add (mpfp_t rop, mpfp_t op1, mpfp_t op2)

{
   mpz_add (rop, op1, op2);
   mpz_mod (rop, rop, __mpfp_p);
}

/*****************************************************************************/

void mpfp_sub (mpfp_t rop, mpfp_t op1, mpfp_t op2)

{
   mpz_sub (rop, op1, op2);
   mpz_mod (rop, rop, __mpfp_p);
}

/*****************************************************************************/

void mpfp_sub_ui (mpfp_t rop, mpfp_t op1, unsigned long int op2)

{
   mpz_sub_ui (rop, op1, op2);
   mpz_mod (rop, rop, __mpfp_p);
}

/*****************************************************************************/

void mpfp_mul (mpfp_t rop, mpfp_t op1, mpfp_t op2)

{
   mpz_mul (rop, op1, op2);
   mpz_mod (rop, rop, __mpfp_p);
   __mpfp_mul++;
}

/*****************************************************************************/

void mpfp_mul_ui (mpfp_t rop, mpfp_t op1, unsigned long int op2)

{
   mpz_mul_ui (rop, op1, op2);
   mpz_mod (rop, rop, __mpfp_p);
   __mpfp_mul++;
}

/*****************************************************************************/

void mpfp_sqr (mpfp_t rop, mpfp_t op)

{
   mpz_pow_ui (rop, op, 2ul);
   mpz_mod (rop, rop, __mpfp_p);
   __mpfp_sqr++;
}

/*****************************************************************************/

void mpfp_add_mul (mpfp_t rop, mpfp_t op1, mpfp_t op2)
   /* sets rop to rop + op1 * op2 */

{
   mpfp_t tmp;

   mpfp_init (tmp);
   mpfp_mul (tmp, op1, op2);
   mpfp_add (rop, rop, tmp);
   mpfp_clear (tmp);
}

/*****************************************************************************/

void mpfp_neg (mpfp_t rop, mpfp_t op)

{
   if (mpz_cmp_ui (op, 0ul) == 0)
      mpz_set (rop, op);
   else
      mpz_sub (rop, __mpfp_p, op);
}

/*****************************************************************************/

void mpfp_inv (mpfp_t rop, mpfp_t op)

{
   mpz_invert (rop, op, __mpfp_p);
   __mpfp_inv++;
}

/*****************************************************************************/
/*****************************************************************************/
