/*

params.h - header file for params.c

Copyright (C) 2009, 2010, 2015, 2021 Andreas Enge

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

#ifndef __PARAMS_H
#define __PARAMS_H

#include <pari/pari.h>
#include <string.h>
#include <getopt.h>
#include "cm.h"

#if defined (__cplusplus)
extern "C" {
#endif

extern void evaluate_parameters (int argc, char* argv [], int_cl_t *d,
   char *invariant, bool *verbose);
extern void evaluate_parameters_ecpp (int argc, char* argv [],
   mpz_ptr n, bool *output, char** filename, bool *verbose, bool *debug,
   bool *check);

#if defined (__cplusplus)
}
#endif
#endif /* ifndef __PARAMS_H */
