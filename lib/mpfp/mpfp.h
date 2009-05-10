#ifndef __MPFP_H
#define __MPFP_H

#include <stdio.h>
#include "gmp.h"

#define MPFPX_MAX_DEG 40
#define MPFPX_INFTY   MPFPX_MAX_DEG

#define MPFPX_NAIVE     0
#define MPFPX_KARATSUBA 1

typedef mpz_t mpfp_t;

#if defined (__cplusplus)
extern "C" {
#endif

extern void mpfp_type_init (mpz_t p);
   /* sets the prime modulus and initialises the "class variables" */

extern void mpfp_init (mpfp_t rop);
extern void mpfp_clear (mpfp_t rop);

extern void mpfp_out (mpfp_t op);

extern void mpfp_set (mpfp_t rop, mpfp_t op);
extern void mpfp_set_z (mpfp_t rop, mpz_t op);
extern void mpfp_set_ui (mpfp_t rop, long unsigned int op);
extern void mpfp_set_str (mpfp_t rop, char *op, int base);

extern void mpfp_init_set (mpfp_t rop, mpfp_t op);
extern void mpfp_init_set_z (mpfp_t rop, mpz_t op);
extern void mpfp_init_set_ui (mpfp_t rop, unsigned long int op);

extern void mpfp_get_z (mpz_t rop, mpfp_t op);

extern int mpfp_cmp_ui (mpfp_t op1, unsigned long int op2);
extern int mpfp_cmp_si (mpfp_t op1, long int op2);

extern void mpfp_add (mpfp_t rop, mpfp_t op1, mpfp_t op2);
extern void mpfp_sub (mpfp_t rop, mpfp_t op1, mpfp_t op2);
extern void mpfp_sub_ui (mpfp_t rop, mpfp_t op1, unsigned long int op2);
extern void mpfp_mul (mpfp_t rop, mpfp_t op1, mpfp_t op2);
extern void mpfp_mul_ui (mpfp_t rop, mpfp_t op1, unsigned long int op2);
extern void mpfp_sqr (mpfp_t rop, mpfp_t op);
extern void mpfp_add_mul (mpfp_t rop, mpfp_t op1, mpfp_t op2);
   /* sets rop to rop + op1 * op2 */
extern void mpfp_neg (mpfp_t rop, mpfp_t op);
extern void mpfp_inv (mpfp_t rop, mpfp_t op);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __MPFP_H */
