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
