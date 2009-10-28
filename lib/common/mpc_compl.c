/*

mpc_compl.c - functions missing in mpc

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

#include "cm_common-impl.h"

/*****************************************************************************/

void mpc_pow_ui (mpc_t rop, mpc_t op1, unsigned long int op2)

{
   mpc_t tmp_op1;

   mpc_init2 (tmp_op1, mpc_get_prec (op1));
   mpc_set (tmp_op1, op1, MPC_RNDNN);

   mpc_set_ui_ui (rop, 1, 0, MPC_RNDNN);
   while (op2 > 0)
   {
      if (op2 % 2 != 0)
      {
         mpc_mul (rop, rop, tmp_op1, MPC_RNDNN);
         op2--;
      }
      mpc_sqr (tmp_op1, tmp_op1, MPC_RNDNN);
      op2 >>= 1;
   }

   mpc_clear (tmp_op1);
}

/*****************************************************************************/
/*****************************************************************************/
