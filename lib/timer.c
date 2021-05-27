/*

timer.c - helper functions for measuring elapsed cpu time

Copyright (C) 2009, 2010, 2021 Andreas Enge

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

#include "cm-impl.h"

/*****************************************************************************/

void cm_timer_start (cm_timer t)

{
   t->elapsed = 0;
   t->time_old = clock ();
}

/*****************************************************************************/

void cm_timer_reset (cm_timer t)

{
   t->elapsed = 0;
}

/*****************************************************************************/

void cm_timer_continue (cm_timer t)

{
   t->time_old = clock ();
}

/*****************************************************************************/

void cm_timer_stop (cm_timer t)

{
   clock_t time_new;
   time_new = clock ();
   t->elapsed += ((double) (time_new - t->time_old)) / CLOCKS_PER_SEC;
}

/*****************************************************************************/

double cm_timer_get (cm_timer t)

{
   return t->elapsed;
}

/*****************************************************************************/
