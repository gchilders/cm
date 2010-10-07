/*

mpfp.h - header file for the mpfp library

Copyright (C) 2009, 2010 Andreas Enge

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

#ifndef __MPFP_H
#define __MPFP_H

#include <stdio.h>
#include "gmp.h"

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
