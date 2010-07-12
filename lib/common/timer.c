/*

timer.c - helper functions for measuring elapsed cpu time

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

#include "cm_common-impl.h"

/*****************************************************************************/

void cm_timer_start (cm_timer clock)

{
   clock->elapsed = 0;
   times (&(clock->time_old));
}

/*****************************************************************************/

void cm_timer_stop (cm_timer clock)

{
   struct tms time_new;
   times (&time_new);
   clock->elapsed =
     ((double) (time_new.tms_utime - clock->time_old.tms_utime)) / 100;
}

/*****************************************************************************/

double cm_timer_get (cm_timer clock)

{
   return clock->elapsed;
}

/*****************************************************************************/
