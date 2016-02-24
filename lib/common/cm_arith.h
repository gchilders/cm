/*

cm_arith.h - macros for arithmetic kernel

Copyright (C) 2015, 2016 Andreas Enge

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

#ifndef __CM_ARITH_H
#define __CM_ARITH_H

#include <mpc.h>

#define ctype              mpc_t
#define cptr               mpc_ptr
#define csrcptr            mpc_srcptr
#define cinit(z,n)         mpc_init2((z),(n))
#define cclear(z)          mpc_clear((z))
#define cset_prec(z,n)     mpc_set_prec((z),(n))
#define cget_prec(z)       mpc_get_prec((z))
#define crealref(z)        mpc_realref((z))
#define cimagref(z)        mpc_imagref((z))
#define cset(z,x)          mpc_set((z),(x),MPC_RNDNN)
#define cset_ui(z,x)       mpc_set_ui((z),(x),MPC_RNDNN)
#define cset_si(z,x)       mpc_set_si((z),(x),MPC_RNDNN)
#define cset_ui_ui(z,x,y)  mpc_set_ui_ui((z),(x),(y),MPC_RNDNN)
#define cneg(z,x)          mpc_neg((z),(x),MPC_RNDNN)
#define cconj(z,x)         mpc_conj((z),(x),MPC_RNDNN)
#define cadd(z,x,y)        mpc_add((z),(x),(y),MPC_RNDNN)
#define cadd_ui(z,x,y)     mpc_add_ui((z),(x),(y),MPC_RNDNN)
#define cadd_si(z,x,y)     mpc_add_si((z),(x),(y),MPC_RNDNN)
#define csub(z,x,y)        mpc_sub((z),(x),(y),MPC_RNDNN)
#define cmul(z,x,y)        mpc_mul((z),(x),(y),MPC_RNDNN)
#define cmul_ui(z,x,y)     mpc_mul_ui((z),(x),(y),MPC_RNDNN)
#define cmul_si(z,x,y)     mpc_mul_si((z),(x),(y),MPC_RNDNN)
#define cmul_fr(z,x,y)     mpc_mul_fr((z),(x),(y),MPC_RNDNN)
#define csqr(z,x)          mpc_sqr((z),(x),MPC_RNDNN)
#define cdiv(z,x,y)        mpc_div((z),(x),(y),MPC_RNDNN)
#define cdiv_ui(z,x,y)     mpc_div_ui((z),(x),(y),MPC_RNDNN)
#define cui_div(z,x,y)     mpc_ui_div((z),(x),(y),MPC_RNDNN)
#define csqrt(z,x)         mpc_sqrt((z),(x),MPC_RNDNN)
#define cnorm(z,x)         mpc_norm((z),(x),MPC_RNDNN)
#define cpow_ui(z,x,y)     mpc_pow_ui((z),(x),(y),MPC_RNDNN)
#define cexp(z,x)          mpc_exp((z),(x),MPC_RNDNN)
#define cinp_str(z,x,y,t)  mpc_inp_str((z),(x),(y),(t),MPC_RNDNN)
#define cout_str(z,x,y,t)  mpc_out_str((z),(x),(y),(t),MPC_RNDNN)

#define ftype              mpfr_t
#define fptr               mpfr_ptr
#define fsrcptr            mpfr_srcptr
#define fprec_t            mpfr_prec_t
#define finit(z,n)         mpfr_init2((z),(n))
#define fclear(z)          mpfr_clear((z))
#define ffree_cache()      mpfr_free_cache()
#define fset_prec(z,n)     mpfr_set_prec((z),(n))
#define fget_prec(z)       mpfr_get_prec((z))
#define fget_exp(z)        mpfr_get_exp((z))
#define fget_z_exp(z,x)    mpfr_get_z_exp((z),(x))
#define fget_si(z)         mpfr_get_si((z),MPFR_RNDN)
#define fget_emin()        mpfr_get_emin()
#define fzero_p(z)         mpfr_zero_p((z))
#define fsgn(z)            mpfr_sgn((z))
#define fcmp_d(z,x)        mpfr_cmp_d((z),(x))
#define fconst_pi(z)       mpfr_const_pi((z),MPFR_RNDN)
#define fset(z,x)          mpfr_set((z),(x),MPFR_RNDN)
#define fset_ui(z,x)       mpfr_set_ui((z),(x),MPFR_RNDN)
#define fset_si(z,x)       mpfr_set_si((z),(x),MPFR_RNDN)
#define fset_z(z,x)        mpfr_set_z((z),(x),MPFR_RNDN)
#define fneg(z,x)          mpfr_neg((z),(x),MPFR_RNDN)
#define fadd(z,x,y)        mpfr_add((z),(x),(y),MPFR_RNDN)
#define fsub(z,x,y)        mpfr_sub((z),(x),(y),MPFR_RNDN)
#define fmul(z,x,y)        mpfr_mul((z),(x),(y),MPFR_RNDN)
#define fmul_si(z,x,y)     mpfr_mul_si((z),(x),(y),MPFR_RNDN)
#define fmul_2ui(z,x,y)    mpfr_mul_2ui((z),(x),(y),MPFR_RNDN)
#define fsqr(z,x)          mpfr_sqr((z),(x),MPFR_RNDN)
#define fdiv(z,x,y)        mpfr_div((z),(x),(y),MPFR_RNDN)
#define fdiv_ui(z,x,y)     mpfr_div_ui((z),(x),(y),MPFR_RNDN)
#define fdiv_2ui(z,x,y)    mpfr_div_2ui((z),(x),(y),MPFR_RNDN)
#define fsqrt(z,x)         mpfr_sqrt((z),(x),MPFR_RNDN)
#define fsqrt_ui(z,x)      mpfr_sqrt_ui((z),(x),MPFR_RNDN)
#define fpow_ui(z,x,y)     mpfr_pow_ui((z),(x),(y),MPFR_RNDN)
#define fexp(z,x)          mpfr_exp((z),(x),MPFR_RNDN)
#define fsin_cos(z,y, x)   mpfr_sin_cos((z),(y),(x),MPFR_RNDN)
#define fround(z,x)        mpfr_round((z),(x))
#define fout_str(z,x,y,t)  mpfr_out_str((z),(x),(y),(t),MPFR_RNDN)

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_ARITH_H */
