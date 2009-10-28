/*

mpfpx.h - header file for the mpfpx library

Copyright (C) 2009 Andreas Enge

This file is part of CM.

CM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the license, or (at your
option) any later version.

CM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with CM; see the file COPYING. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#ifndef __MPFPX_H
#define __MPFPX_H

#include <stdio.h>
#include "mpfp.h"

#define MPFPX_MAX_DEG 40
#define MPFPX_INFTY   MPFPX_MAX_DEG

#define MPFPX_NAIVE     0
#define MPFPX_KARATSUBA 1

typedef struct
{
   int    deg;
   mpfp_t coeff [MPFPX_MAX_DEG + 1];
}
__mpfpx_t;
typedef __mpfpx_t mpfpx_t [1];


#if defined (__cplusplus)
extern "C" {
#endif

extern void mpfpx_type_init (mpz_t p, int arithmetic);
   /* sets the prime modulus and initialises the "class variables" */

extern void mpfpx_init (mpfpx_t rop);
extern void mpfpx_clear (mpfpx_t rop);
extern void mpfpx_normalise (mpfpx_t rop);
   /* deletes leading zeroes */
extern void mpfpx_monicise (mpfpx_t rop, mpfpx_t op, mpfp_t lc_inv);
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op.                                                */

extern void mpfpx_set (mpfpx_t rop, mpfpx_t op);
extern void mpfpx_set_fp (mpfpx_t rop, mpfp_t op);
extern void mpfpx_set_ui (mpfpx_t rop, unsigned long int op);
extern void mpfpx_set_z_array (mpfpx_t rop, mpz_t *op, int size);
extern void mpfpx_set_ui_array (mpfpx_t rop, const unsigned long int *op,
   int size);
extern void mpfpx_set_str_array (mpfpx_t rop, char **op, int base, int size);

extern void mpfpx_init_set (mpfpx_t rop, mpfpx_t op);
extern void mpfpx_init_set_fp (mpfpx_t rop, mpfp_t op);
extern void mpfpx_init_set_ui (mpfpx_t rop, unsigned long op);
extern void mpfpx_init_set_ui_array (mpfpx_t rop, const unsigned long int *op,
   int size);
extern void mpfpx_init_set_str_array (mpfpx_t rop, char **op, int base,
   int size);

extern void mpfpx_out (mpfpx_t op);

/* The following functions change the degree and coefficients of the       */
/* polynomial; it is in the programmer's responsability to ensure that the */
/* result is coherent.                                                     */
extern void mpfpx_deg_set (mpfpx_t rop, int op);
extern void mpfpx_coeff_set_z (mpfpx_t rop, unsigned int op1, mpz_t op2);
extern void mpfpx_coeff_set_ui (mpfpx_t rop, unsigned int op1,
   unsigned long int op2);
extern void mpfpx_coeff_get_z (mpz_t rop, mpfpx_t op1, unsigned int op2);

extern void mpfpx_add (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2);
extern void mpfpx_sub (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2);
extern void mpfpx_sub_ui (mpfpx_t rop, mpfpx_t op1, unsigned long int op2);
extern void mpfpx_mul (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2);
extern void mpfpx_mul_ui (mpfpx_t rop, mpfpx_t op1, unsigned long int op2);
extern void mpfpx_mul_low (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp);
   /* sets the degree of rop to that of op1*op2 and sets the coefficients */
   /* of degree at most exp of the product                                */
extern void mpfpx_mul_high (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2, char exp);
   /* ditto for the hight degree coefficients, of which we compute exp+1 */
extern void mpfpx_mul_m (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2);
   /* assumes that op2 is monic */
extern void mpfpx_mul_m_low (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2,
   char exp);
extern void mpfpx_mul_m_high (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2,
   char exp);
extern void mpfpx_mul_mm_high (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2,
   char exp);
   /* assumes that op1 and op2 are monic */
extern void mpfpx_sqr (mpfpx_t rop, mpfpx_t op);
extern void mpfpx_sqr_low (mpfpx_t rop, mpfpx_t op, char exp);
extern void mpfpx_sqr_high (mpfpx_t rop, mpfpx_t op, char exp);
extern void mpfpx_sqr_m_high (mpfpx_t rop, mpfpx_t op, char exp);
extern void mpfpx_div_qr (mpfpx_t q, mpfpx_t r, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */
extern void mpfpx_div_qr_newton (mpfpx_t q, mpfpx_t r, mpfpx_t op1,
   mpfpx_t op2);
   /* ditto, but using Newton iterations; in this case, op2 must be monic! */
extern void mpfpx_div_q (mpfpx_t q, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */
extern void mpfpx_div_r (mpfpx_t q, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */
extern void mpfpx_divexact_jebelean (mpfpx_t q, mpfpx_t op1, mpfpx_t op2);
extern void mpfpx_divexact (mpfpx_t q, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */
extern void mpfpx_neg (mpfpx_t rop, mpfpx_t op);

extern void mpfpx_gcd_ext (mpfpx_t d, mpfpx_t u, mpfpx_t v, mpfpx_t a,
   mpfpx_t b, int steps, mpfp_t lc_inv);
   /* computes at most "steps" steps of the extended Euclidian algorithm and  */
   /* returns the result in d. If "steps" equals CM_MPFPX_INFTY, then         */
   /* d = gcd (a, b).                                                         */
   /* If u and v do not equal NULL, they are computed to yield d = u a + v b. */
   /* If lc_inv is not NULL, it must contain the inverse of the leading       */
   /* coefficient of b.                                                       */
extern void mpfpx_invert (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* computes the inverse of op1 modulo op2                            */
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __MPFPX_H */
