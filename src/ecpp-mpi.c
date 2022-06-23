/*

ecpp-cm.c - executable for elliptic curve primality proofs with MPI

Copyright (C) 2021, 2022 Andreas Enge

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

#include <mpi.h>
#include "params.h"

int main (int argc, char* argv [])
{
   int rank;

   cm_pari_init ();
   cm_mpi_init (true);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if (rank == 0) {
      mpz_t n;
      bool print, verbose, debug, trust, check, onlys1;
      char *filename;
      mpz_init (n);
      evaluate_parameters_ecpp (argc, argv, n, &print, &filename,
         &verbose, &debug, &trust, &check, &onlys1);
      cm_ecpp (n, CM_MODPOLDIR, print, filename, trust, check, verbose,
         debug, onlys1);
      mpz_clear (n);
   }
   cm_mpi_clear ();
   cm_pari_clear ();

   return 0;
}
