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

void cm_timer_reset (cm_timer_t t)

{
   t->elapsed = 0;
   t->wc_elapsed = 0;
}

/*****************************************************************************/

void cm_timer_continue (cm_timer_t t)

{
   t->time_old = clock ();
   gettimeofday (t->wc_time_old, NULL);
}

/*****************************************************************************/

void cm_timer_start (cm_timer_t t)

{
   cm_timer_reset (t);
   cm_timer_continue (t);
}

/*****************************************************************************/

void cm_timer_stop (cm_timer_t t)

{
   clock_t time_new;
   struct timeval wc_time_new [1];

   time_new = clock ();
   t->elapsed += ((double) (time_new - t->time_old)) / CLOCKS_PER_SEC;
   gettimeofday (wc_time_new, NULL);
   t->wc_elapsed +=   (double) wc_time_new->tv_sec
                    - (double) t->wc_time_old->tv_sec
                    + (  (double) wc_time_new->tv_usec
                       - (double) t->wc_time_old->tv_usec) / 1e6;
}

/*****************************************************************************/

double cm_timer_get (cm_timer_t t)

{
   return t->elapsed;
}

/*****************************************************************************/

double cm_timer_wc_get (cm_timer_t t)

{
   return t->wc_elapsed;
}

/*****************************************************************************/
/*****************************************************************************/

void cm_stat_init (cm_stat_t stat)
{
   unsigned long int i;

   for (i = 0; i < sizeof (stat->counter) / sizeof (stat->counter [0]);
        i++)
      stat->counter [i] = 0;
   for (i = 0; i < sizeof (stat->timer) / sizeof (stat->timer [0]); i++)
      cm_timer_reset (stat->timer [i]);
}

/*****************************************************************************/
/*****************************************************************************/
