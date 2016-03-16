/*

mpfpx.h - header file for the mpfpx library

Copyright (C) 2009, 2010 Andreas Enge

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

#ifndef __MPFPX_H
#define __MPFPX_H

#include <stdio.h>
#include "gmp.h"

#define MPFPX_MAX_DEG 40

#define MPFPX_NAIVE     0
#define MPFPX_KARATSUBA 1

typedef mpz_t mpfp_t;

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

extern void mpfp_init (mpfp_t rop);
extern void mpfp_clear (mpfp_t rop);
extern void mpfp_set_ui (mpfp_t rop, long unsigned int op);


extern void mpfpx_type_init (mpz_t p, int arithmetic);
   /* sets the prime modulus and initialises the "class variables" */

extern void mpfpx_init (mpfpx_t rop);
extern void mpfpx_clear (mpfpx_t rop);

extern void mpfpx_set_z_array (mpfpx_t rop, mpz_t *op, int size);
extern void mpfpx_set_ui_array (mpfpx_t rop, const unsigned long int *op,
   int size);

extern void mpfpx_out (mpfpx_t op);

/* The following functions change the degree and coefficients of the       */
/* polynomial; it is in the programmer's responsability to ensure that the */
/* result is coherent.                                                     */
extern void mpfpx_coeff_get_z (mpz_t rop, mpfpx_t op1, unsigned int op2);

extern void mpfpx_sub_ui (mpfpx_t rop, mpfpx_t op1, unsigned long int op2);
extern void mpfpx_mul (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2);
extern void mpfpx_mul_ui (mpfpx_t rop, mpfpx_t op1, unsigned long int op2);
extern void mpfpx_sqr (mpfpx_t rop, mpfpx_t op);
extern void mpfpx_div_r (mpfpx_t q, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */

extern void mpfpx_invert (mpfpx_t rop, mpfpx_t op1, mpfpx_t op2,
   mpfp_t lc_inv);
   /* computes the inverse of op1 modulo op2                            */
   /* If lc_inv is not NULL, it must contain the inverse of the leading */
   /* coefficient of op2.                                               */

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __MPFPX_H */
