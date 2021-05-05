/*

cm_common-impl.h - header file for internal use of the cm_common library

Copyright (C) 2009, 2012, 2015, 2021 Andreas Enge

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

#ifndef __CM_COMMON_IMPL_H
#define __CM_COMMON_IMPL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "cm_common.h"

#define CM_QDEV_ETA   1
#define CM_QDEV_ATKIN 2


#if defined (__cplusplus)
extern "C" {
#endif

/* functions for computing q expansions of modular functions and addition
   chains */
extern void cm_qdev_init (cm_qdev_t *f, fprec_t prec);
extern void cm_qdev_clear (cm_qdev_t *f);
extern void cm_qdev_eval (ctype rop, cm_qdev_t f, ctype q1);
extern void cm_qdev_eval_fr (ftype rop, cm_qdev_t f, ftype q1);

/* functions for evaluating modular functions */
extern void cm_modular_eta_series_fr (cm_modular_t m, ftype rop, ftype q_24);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __CM_COMMON_IMPL_H */