#include "cm_common-impl.h"

/*****************************************************************************/

void cm_timer_start (cm_timer clock)

{
   clock->elapsed = 0;
   times (&(clock->time_old));
}

/*****************************************************************************/

void cm_timer_run (cm_timer clock)

{
   times (&(clock->time_old));
}

/*****************************************************************************/

void cm_timer_stop (cm_timer clock)

{
   times (&(clock->time_new));
   clock->elapsed +=
     ((double) (clock->time_new.tms_utime - clock->time_old.tms_utime)) / 100;
}

/*****************************************************************************/

double cm_timer_get (cm_timer clock)

{
   return clock->elapsed;
}

/*****************************************************************************/
